#include <iostream>
#include <deque>
#include <mpi.h>
#include "option.h"
#include "logger.h"
#include "basic_types.h"
#include "mpi_interface.h"
#include "file_reader.h"
#include "exchange_buffer.h"
#include "kdtree_balanced.h"
#include "nbr_search.h"
#include "fof.h"
#include "open.h"

using namespace std;

void for_each_readshift(const char redshift[],
			const char particle_dir[],
			const bool shift,
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			const float buffer_factor,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer);


int remove_extra_particles(Particles* const particles, const float buffer_factor);


int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  MpiInterface* const mpi= new MpiInterface();

  //
  // MPI setup
  //

  if(argc == 1 && mpi->index() == 0) {
    cout << "clumping_tree [options]\n"
	 << "\t-pm_redshift <z>; p3m file <z>xv<i>.dat\n"
	 << "\t-pdir <particle directory>; -pdir particle_data\n"
	 << "\t-allocate <MB for particles>\n"
	 << "\t-boxsize <boxsize>\n"
	 << "\t-buffer_factor 0.05\n"
         << "\t-l 0.2; FOF linking length"
         << "\t-nc <nc1,nc2,...> list of number of mesh";
    
    return 0;
  }

  //
  // Begin Logging
  //
  Logger logger;


  //
  // Commandline Options
  //
  Option op(argc, argv);
  op.set_default("-pm_redshift", "redshifts.txt");
  op.set_default("-pdir", "particle_data");
  op.set_default("-hdir", "halo_data");
  op.set_default("-allocate", "512"); // 12M particles ~ 512MB
  op.set_default("-boxsize", "-1");
  op.set_default("-buffer_factor", "0.05");
  op.set_default("-l", "0.2");
  op.set_default("-ncf", "0"); // finemesh clumping
  op.set_default("-shift", "0");

  const float boxsize= op.get_float("-boxsize");
  assert(boxsize > 0.0f);
  const float buffer_factor= op.get_float("-buffer_factor");
  assert(buffer_factor > 0.0f);
  //const int nc= op.get_int("-nc");
  //const bool full_box= 0; //op.get_int("-full");
  const float ll= 2.0f * op.get_float("-l");
  const bool do_shift= op.get_int("-shift");

  logger << "boxsize " << boxsize << "\n";
  logger << "buffer_factor " << buffer_factor << "\n";
  logger << "linking_length(internal) " << ll << "\n";
  
  // open redshift list
  string redshift_filename= op["-pm_redshift"];
  deque<string> redshifts;
  if(!redshift_filename.empty()) {
    FILE* fp= fopen(redshift_filename.c_str(), "r");
    assert(fp);
    
    char buf[128], z[32];
    while(fgets(buf, 128, fp)) {
      if(sscanf(buf, "%s", z) == 1) {
	redshifts.push_back(string(z));
      }
    }

    fclose(fp);
  }
  else {
    mpi->abort("No redshift file");
  }

  //
  // Allocate Memory
  //
  const size_t mbyte_alloc= op.get_int("-allocate");
  const size_t np_allocate= mbyte_alloc*1024*1024/sizeof(Particle);
  Particles* particles= new Particles();
  particles->allocate(np_allocate);
  particles->boxsize= boxsize;

  KDTree* tree= new KDTreeSimple();
  tree->allocate(np_allocate, 16);
  set_nbr_search(tree->root(), particles->particle);

  int* const temp= (int*) malloc(sizeof(float)*np_allocate);
  assert(temp);

  ParticleSet<Halo>* halos= new ParticleSet<Halo>();
  halos->allocate(temp, np_allocate/sizeof(Halo));
  //halos->allocate(2000);
  halos->boxsize= boxsize;

  //
  // Make output directory
  //
  if(mpi->rank() == 0) {
    char dirname[128];
    int ret= make_directory("fof_halos/"); assert(ret);
  }


  //
  // For each snapshot
  //
  for(deque<string>::iterator z= redshifts.begin();
      z != redshifts.end(); ++z) {
    for_each_readshift(z->c_str(),        // redshift
		       op["-pdir"].c_str(), // particle dir
		       do_shift,
		       mpi,
		       logger,
		       particles,
		       buffer_factor,
		       ll,
		       tree,
		       temp);

  }

  logger << "done all\n";

  delete tree;
  delete mpi;
  free(temp);

  MPI_Finalize();

  return 0;
}

  /*
  //
  // Write rk
  //
  FILE* fp= fopen("rk1", "w");
  assert(fp);
  fwrite(&count, sizeof(int), 1, fp);
  for(index_t i=0; i<np; ++i) {
    if(0 <= p[i].x[0] && p[i].x[0] < boxsize &&
       0 <= p[i].x[1] && p[i].x[1] < boxsize &&
       0 <= p[i].x[2] && p[i].x[2] < boxsize) {
      fwrite(&(p[i].id), sizeof(int), 1, fp);
      fwrite(&(p[i].rk), sizeof(float), 1, fp);
    }
  }
  fwrite(&count, sizeof(int), 1, fp);
  fclose(fp);
  */


void for_each_readshift(const char redshift[],
			const char particle_dir[],
			const bool do_shift,
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			const float buffer_factor,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer)
{
  const float boxsize= particles->boxsize;
  //
  // Read File
  //
  mpi->print("reading particles\n");
  logger.begin_timer(io);

  logger << "z " << redshift << "\n";
  char filename[256];

  //
  // Read and distribute shift
  //
  float shift[]= {0.0f, 0.0f, 0.0f};
  if(do_shift) {
    if(mpi->rank() == 0) {
      sprintf(filename, "shift/%sshift.txt", redshift);
      float z;
      FILE* fp= fopen(filename, "r");
      if(fp == 0) {
	cerr << "Unable to open shift: " << filename << endl;
	mpi->abort("shift file error");
      }
      int ret= fscanf(fp, "%e %e %e %e", &z, shift, shift+1, shift+2);
      assert(ret == 4);
      fclose(fp);
    }
    mpi->bcast(shift, 3, MPI_FLOAT);
  }
  logger << "shift " 
	 << shift[0] << " " << shift[1] << " " << shift[2] << "\n";

  if(mpi->rank() == 1) { // debug
    cout << "shift(1) " << shift[0] << " " << shift[1] << " " << shift[2] << "\n";
  }
  
  sprintf(filename, "%s/%sxv%d.dat", 
	  particle_dir, redshift, mpi->index());

  read_pm_file2(filename, particles, shift);

  if(particles->np_local <= 0)
    mpi->abort("No particles\n");

  logger << "np_local " << particles->np_local << "\n";
  MPI_Barrier(MPI_COMM_WORLD);


  //
  // Exchange buffer particles
  //
  mpi->print("exchanging buffer particles\n");
  logger.begin_timer(communication);
  exchange_buffer(mpi, particles, buffer_factor);
  logger << "np_with_buffers1 " << particles->np_with_buffers << "\n";

  //
  // Remove Extra particles outside buffer zone
  //
  mpi->print("removing extra particles\n");
  int np_removed= remove_extra_particles(particles, buffer_factor);
  logger << "np_removed " << np_removed << "\n";
  logger << "np_with_buffers " << particles->np_with_buffers << "\n";

  //
  // Construct KDTree
  //
  mpi->print("constructing kdtree\n");
  logger.begin_timer(construction);
  //KDTree* tree= new KDTreeBalanced();
  const index_t nnode= tree->construct(particles, 
				 -buffer_factor*boxsize, 
				 (1.0f+buffer_factor)*boxsize);
  logger << "nnode " << nnode << "\n";


  //
  // FOF
  //
  mpi->print("fof\n");
  logger.begin_timer(fof);

  const index_t ngrp= group_fof(particles, tree, ll);
  logger << "ngrp " << ngrp << "\n";
  const index_t nhalo= remove_small_groups(particles, temp_buffer, ngrp);
  logger << "n_fof_halo " << nhalo << "\n";
  
  sprintf(filename, "fof_halos/%sfof%d.dat", redshift, mpi->index());
  write_fof_halos(filename, particles, temp_buffer, nhalo);

  logger.begin_timer(imbalance);
  MPI_Barrier(MPI_COMM_WORLD);

  logger.end_timer();
    
  if(mpi->index() == 0)
    logger.print_log(cout);

  mpi->print("bye\n\n");
}

int remove_extra_particles(Particles* const particles, const float buffer_factor)
{
  index_t np= particles->np_with_buffers;
  Particle* const p= particles->particle;
  Particle* q= particles->particle;
  
  const float dx= buffer_factor*particles->boxsize;
  const float left= -dx;
  const float right= particles->boxsize+dx;
  
  index_t np_after= 0;

  for(int i=0; i<np; ++i) {
    if(left < p[i].x[0] && p[i].x[0] < right &&
       left < p[i].x[1] && p[i].x[1] < right &&
       left < p[i].x[2] && p[i].x[2] < right) {
      *q= p[i];
      q++;
      np_after++;
    }
  }

  particles->np_local -= (np - np_after);
  particles->np_with_buffers= np_after;

  return (np - np_after);
}

