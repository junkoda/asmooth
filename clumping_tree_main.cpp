#include <cstring>
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
#include "coarse_mesh.h"
#include "halo_file.h"
#include "open.h"

using namespace std;

void for_each_readshift(const char redshift[],
			const char particle_dir[],
			const char halo_dir[],
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			ParticleSet<Halo>* halos,
			const float buffer_factor,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer,
			float* const mesh,
			deque<int>& ncs);

int remove_extra_particles(Particles* const particles, const float buffer_factor);

void check_total_mass(const float mesh[], const int nc, const long np_total, const int mpi_rank);

char* pop_num(char* const p, int* const n)
{
  char* q= strchr(p, ',');
  if(q != 0)
    *q= 0;
  
  *n= atoi(p);

  return q == 0 ? 0 : q+1;
}

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
  op.set_default("-nc", "256");       // nc^3 is number of mesh in this node
  //op.set_default("-full", "1");       // shift halo coordinate
  op.set_default("-ncf", "0"); // finemesh clumping

  const float boxsize= op.get_float("-boxsize");
  assert(boxsize > 0.0f);
  const float buffer_factor= op.get_float("-buffer_factor");
  assert(buffer_factor > 0.0f);
  //const int nc= op.get_int("-nc");
  //const bool full_box= 0; //op.get_int("-full");
  const float ll= 2.0f * op.get_float("-l");

  logger << "boxsize " << boxsize << "\n";
  logger << "buffer_factor " << buffer_factor << "\n";
  logger << "linking_length(internal) " << ll << "\n";
  
  // nc list
  char nc_arg[128];
  char * nc_list= nc_arg;
  strncpy(nc_arg, op["-nc"].c_str(), 127);
  deque<int> ncs;
  while(nc_list != 0) {
    int nc1= 0;
    nc_list= pop_num(nc_list, &nc1);
    ncs.push_back(nc1);
  }
  int nc= 0;
  logger << "nc ";
  for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
    logger <<*p << " ";
    if(*p > nc) nc= *p;
  }
  logger << "\n";
  assert(nc > 0);
  
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

  float* const mesh= (float*) malloc(sizeof(float)*nc*nc*nc*6);
  assert(mesh);

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

    for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
      const int nc= *p; assert(nc > 0);

      sprintf(dirname, "fof/nc%d/", nc);
      ret= make_directory(dirname); assert(ret);

      sprintf(dirname, "so/nc%d/", nc);
      make_directory(dirname); assert(ret);
    }
  }


  //
  // For each snapshot
  //
  for(deque<string>::iterator z= redshifts.begin();
      z != redshifts.end(); ++z) {
    for_each_readshift(z->c_str(),        // redshift
		       op["-pdir"].c_str(), // particle dir
		       op["-hdir"].c_str(), // halo dir
		       mpi,
		       logger,
		       particles,
		       halos,
		       buffer_factor,
		       ll,
		       tree,
		       temp,
		       mesh,
		       ncs);

  }

  logger << "done all\n";

  delete tree;
  delete mpi;
  free(mesh);
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
			const char halo_dir[],
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			ParticleSet<Halo>* halos,
			const float buffer_factor,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer,
			float* const mesh,
			deque<int>& ncs)
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
  logger << "shift " 
	 << shift[0] << " " << shift[1] << " " << shift[2] << "\n";

  if(mpi->rank() == 1) { // debug
    cout << "shift(1) " << shift[0] << " " << shift[1] << " " << shift[2] << "\n";
  }
  
  //string pm_redshift= op["-pm_redshift"];
  //if(!pm_redshift.empty()) {    
  sprintf(filename, "%s/%sxv%d.dat", 
	  particle_dir, redshift, mpi->index());
    	    //op["-pdir"].c_str(), pm_redshift.c_str(), mpi->rank());
    //	    op["-pdir"].c_str(), pm_redshift.c_str(), 0); // debug

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

  //
  // Neighbor Search
  //
  
  //const int k_nbr= 32;
  mpi->print("neighbor search\n");
  logger.begin_timer(neighbors);
  KNeighbors knbrs;
  //set_nbr_search(tree->root(), particles->particle);
  
  Particle* const p= particles->particle;
  const index_t np= particles->np_with_buffers;
  //  int count= 0;

#ifdef _OPENMP
#pragma omp parallel for private(knbrs)
#endif
  for(index_t i=0; i<np; ++i) {
    neighbor_search(p+i, &knbrs);
  }


  //
  // Mesh (FOF)
  //
  for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
    const int nc= *p; assert(nc > 0);

    logger.begin_timer(mesh_assign);
    clear_mesh(mesh, nc);
    assign_on_mesh(particles, mesh, nc); // local
    check_total_mass(mesh, nc, particles->np_total, mpi->rank());

    logger.begin_timer(output);
    sprintf(filename, "fof/nc%d/%s", nc, redshift);
    write_mesh_separate(filename, mpi->index(), mesh, nc, boxsize);
  }
  //
  //sprintf(filename, "fof/%smesh%d.dat", redshift, mpi->index());
  //write_mesh(filename, mesh, nc, boxsize);

  //logger.begin_timer(imbalance);
  MPI_Barrier(MPI_COMM_WORLD);

  //
  // Spherical Overdensity Halo
  //
  logger.begin_timer(so_halo);
  sprintf(filename, "%s/%shalo%d.dat", 
	  halo_dir, redshift, mpi->index());
  read_and_exchage_halo(filename, mpi, halos, buffer_factor, shift);
  Halo const * const h= halos->particle;
  const index_t n_sohalo= halos->np_with_buffers;

  // Mark Halo Particles
  const float left= -boxsize * buffer_factor;
  const float right= boxsize + boxsize * buffer_factor;
  long np_halo= 0;
  double m_halo= 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:np_halo, m_halo) shared(h)
#endif
  for(index_t i=0; i<n_sohalo; ++i) {
    //cerr << "rhalo " << h[i].r << endl;
    if(left < h[i].x[0] && h[i].x[0] < right &&
       left < h[i].x[1] && h[i].x[1] < right &&
       left < h[i].x[2] && h[i].x[2] < right) {
      int nph= mark_halo_particles(tree->root(), h[i].x, h[i].r);
      np_halo += nph;
      m_halo += h[i].mass/8.0f;
      /*
      if(mpi->index() == 0) // debug!!!
	printf("hdebug %e %e %e %e %d\n", h[i].x[0], h[i].x[1], h[i].x[2], 
	       h[i].mass/8.0f, nph);
      */
    }
  }

  logger << "n_so_halo " << halos->np_local << " " << halos->np_with_buffers << "\n"; 

  logger << "so_halo_local " << np_halo << " " << m_halo << "\n";

  logger.begin_timer(imbalance);
  MPI_Barrier(MPI_COMM_WORLD);

  // Check total halo mass
  logger.begin_timer(communication);

  long np_halo_total= 0;
  double m_halo_total= 0.0f;
  MPI_Reduce(&np_halo, &np_halo_total, 1, MPI_LONG, 
	     MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&m_halo, &m_halo_total, 1, MPI_DOUBLE, 
	     MPI_SUM, 0, MPI_COMM_WORLD);

  if(mpi->rank() == 0) {
    cout << "check halo_mass " << redshift << " " 
	 << np_halo_total << " " << m_halo_total << "\n";
    if(m_halo_total > 100.0) {
      assert(0.80 < np_halo_total/m_halo_total && 
	     np_halo_total/m_halo_total < 1.20);
    }
  }

  // debug !
  /*
  np_halo= 0;
  for(index_t i=0; i<np; ++i)
    np_halo += (particles->particle[i].igrp == -1);

  logger << "n_so_halo_local " << halos->np_local << "\n";
  logger << "np_so_halo " << np_halo << " (with buffers) \n";
  */

  //  
  // Mesh Assignment (SO)
  //

  for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
    logger.begin_timer(mesh_assign);

    const int nc= *p;
    clear_mesh(mesh, nc);
    assign_on_mesh(particles, mesh, nc, -1);

    check_total_mass(mesh, nc, particles->np_total, mpi->rank());

    logger.begin_timer(output);
    sprintf(filename, "so/nc%d/%s", nc, redshift);
    write_mesh_separate(filename, mpi->index(), mesh, nc, boxsize);
  }
  //assign_on_global_mesh(particles, mpi, mesh, nc, -1);
  //write_mesh("sog/mesh.dat", mesh, nc, boxsize);




  //
  // Debug! (test using finemesh)
  //
  /*
  logger.begin_timer(test);
  const int ncf= op.get_int("-ncf");
  if(ncf > 0) {
    const int nfmesh= ncf*ncf*ncf;
    float* finemesh= (float*) malloc(sizeof(float)*nfmesh);
    assign_on_density_mesh(particles, finemesh, ncf, -1);
    double nf2= 0.0;
    float vf_inv= 1.0*nfmesh/(boxsize*boxsize*boxsize);
    for(int i=0; i<nfmesh; ++i) {
      float n= finemesh[i]*vf_inv;
      nf2 += n*n;
    }
    
    logger << "nf2 " << ncf << " " << nf2/nfmesh*64.0 << endl;
  }
  */
  // nbar= 1/8.0; 1/nbar^2= 64.0
  
  //
  // Done
  //

  logger.begin_timer(imbalance);
  MPI_Barrier(MPI_COMM_WORLD);

  logger.end_timer();
    
  if(mpi->index() == 0)
    logger.print_log(cout);

  mpi->print("bye\n");
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

void check_total_mass(const float mesh[], const int nc, const long np_total, const int mpi_rank)
{
  // Check total mass
  double m_sum[2], m_global[2];
  m_sum[0]= mesh_total(mesh, nc);
  m_sum[1]= mesh_total(mesh+nc*nc*nc, nc);

  //logger.begin_timer(imbalance);
  MPI_Barrier(MPI_COMM_WORLD);

  //logger.begin_timer(communication);
  MPI_Reduce(m_sum, m_global, 2, MPI_DOUBLE, 
	     MPI_SUM, 0, MPI_COMM_WORLD);
  double m_total= m_global[0] + m_global[1];

  if(mpi_rank == 0) {
    printf("check total_mass %4d %15.1f %15.1f %ld\n",
	   nc, m_global[0], m_total, np_total);
    
    assert(0.999 < m_total/np_total && 
	           m_total/np_total < 1.001);
  }
}
