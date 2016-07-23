#include <cstring>
#include <iostream>
#include <deque>
#include <mpi.h>
#include "option.h"
#include "logger.h"
#include "basic_types.h"
#include "mpi_interface.h"
#include "file_reader2.h"
#include "exchange_buffer.h"
#include "kdtree_balanced.h"
#include "nbr_search.h"
#include "fof.h"
#include "coarse_mesh.h"
#include "halo_file.h"
#include "open.h"

using namespace std;

// Semi-global variables
static float boxsize;
static bool read_xv_file= false;

static void for_each_redshift(const char redshift[],
			const char filebase[],
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			ParticleSet<Halo>* halos,
			const float buffer_factor,
			const int nc_node_dim, const int mesh_scale,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer,
			float* const mesh,
			deque<int>& ncs);


static int remove_extra_particles(Particles* const particles, const float buffer_factor);

static void check_total_mass(const float mesh[], const int nc, const long np_total, const int mpi_rank);

static char* pop_num(char* const p, int* const n)
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
	 << "\t-pm_redshift <redshifts.txt>; file with list of redshifts\n"
	 << "\t-node_dir <node directory >; -node_dir results/node\n"
	 << "\t-allocate <MB for particles>\n"
	 << "\t-buffer_factor 0.05\n"
	 << "\t-xv 912; -xv <local boxsize> (for xv particle format)\n"
         << "\t-nc_node_dim 32 (for zip particle format)\n"
	 << "\t-mesh_scale 4   (for zip particle format)\n"
         << "\t-l 0.2; FOF linking length\n"
         << "\t-omegam 0.308; (for Bryan-Normal over density\n"
         << "\t-nc <nc1,nc2,...> list of number of mesh\n";
    
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
  op.set_default("-node_dir", "results/node");
  //op.set_default("-hdir", "halo_data");
  op.set_default("-allocate", "512"); // 12M particles ~ 512MB
  op.set_default("-xv", "-1");
  op.set_default("-buffer_factor", "0.05");
  op.set_default("-nc_node_dim", "0");
  op.set_default("-l", "0.2");
  op.set_default("-nc", "256");       // nc^3 is number of mesh in this node
  op.set_default("-mesh_scale", "4");
  op.set_default("-omegam", "0");

  // If you read particles from xv file, give boxsize with -xv boxsize
  // boxsize is in internal unit; 2xnumber of particles per dim
  boxsize= op.get_float("-xv");
  
  const float buffer_factor= op.get_float("-buffer_factor");
  const int nc_node_dim= op.get_int("-nc_node_dim");
  const int mesh_scale= op.get_int("-mesh_scale");

  assert(buffer_factor > 0.0f);
  const float ll= 2.0f * op.get_float("-l");

  logger << "buffer_factor " << buffer_factor << "\n";
  logger << "linking_length(internal) " << ll << "\n";

  if(boxsize <= 0.0f) {
    read_xv_file = false;
    logger << "zip particle format\n";
    logger << "nc_node_dim " << nc_node_dim << "\n";
    logger << "mesh_scale " << mesh_scale << "\n";

    if(nc_node_dim <= 0) {
      mpi->abort("Error: option -nc_node_dim not given");
    }
  }
  else {
    read_xv_file = true;
    logger << "-xv option given; boxsize = " << boxsize << "\n";
  }
  
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
    if(fp == 0) {
      cerr << "Error: Redshift file " << redshift_filename << " does not exist\n";
    }
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

  if(mpi->rank() == 0) {
    cerr << "MPI size = " << mpi->size() << endl;
    cerr << "Allocate " << mbyte_alloc << " Mbytes for particles\n";
    cerr << "np_allocate = " << np_allocate << endl;
  }

  Particles* particles= new Particles();
  particles->allocate(np_allocate);
  particles->omega_m= op.get_float("-omegam");

  KDTree* tree= new KDTreeSimple();

  size_t mem_tree= tree->estimate_mem(np_allocate, 16);
  if(mpi->rank() == 0) {
    cerr << mem_tree << " bytes for KD tree" << endl;
    cerr << mem_tree/(1024*1024) << " Mbytes\n";
  }
  
  tree->allocate(np_allocate, 16);
  set_nbr_search(tree->root(), particles->particle);

  if(mpi->rank() == 0) {
    cerr << sizeof(float)*np_allocate << " bytes for temp" << endl;
    cerr << "or, " << sizeof(float)*np_allocate/(1024*1024) << " Mbytes\n";
  }


  int* const temp= (int*) malloc(sizeof(float)*np_allocate);
  assert(temp);

  if(mpi->rank() == 0) {
    cerr << (sizeof(float)*nc*nc*nc*6) << " bytes for mesh" << endl;
    cerr << "or, " << (sizeof(float)*nc*nc*nc*6/(1024*1024)) << " Mbytes\n";
  }
  
  float* const mesh= (float*) malloc(sizeof(float)*nc*nc*nc*6);
  assert(mesh);

  ParticleSet<Halo>* halos= new ParticleSet<Halo>();
  halos->allocate(temp, np_allocate/sizeof(Halo));

  //
  // Make output directory
  //
  if(mpi->rank() == 0) {
    char dirname[128];
    int ret= make_directory("fof_halos/"); assert(ret);

    for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
      const int nc= *p; assert(nc > 0);
#ifdef FOFON
      sprintf(dirname, "fof/nc%d/", nc);
      ret= make_directory(dirname); assert(ret);
#endif

      sprintf(dirname, "so/nc%d/", nc);
      make_directory(dirname); assert(ret);
    }
  }


  //
  // For each snapshot
  //
  for(deque<string>::iterator z= redshifts.begin();
      z != redshifts.end(); ++z) {
    for_each_redshift(z->c_str(),        // redshift
		       op["-node_dir"].c_str(), // particle dir
		       mpi,
		       logger,
		       particles,
		       halos,
		       buffer_factor,
		       nc_node_dim, mesh_scale,
		       ll,
		       tree,
		       temp,
		       mesh,
		       ncs);

  }

  logger << "boxsize " << particles->boxsize << "\n";
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


void for_each_redshift(const char redshift[],
			const char filebase[],
			MpiInterface const * const mpi,
			Logger& logger,
			Particles* const particles,
			ParticleSet<Halo>* halos,
			const float buffer_factor,
			const int nc_node_dim, const int mesh_scale,
			const float ll,
			KDTree* const tree,
			int* const temp_buffer,
			float* const mesh,
			deque<int>& ncs)
{
  //
  // Read File
  //
  mpi->print("reading particles\n");
  logger.begin_timer(io);

  logger << "z " << redshift << "\n";
  const float z= atof(redshift);
  const float omega_m= particles->omega_m;

  //
  // Read and distribute shift
  //
  char filename[256];
  int shift_file= true;
  float shift[]= {0.0f, 0.0f, 0.0f};

  if(mpi->rank() == 0) {
    sprintf(filename, "shift/%sshift.txt", redshift);
    float z;
    FILE* fp= fopen(filename, "r");
    if(fp == 0) {
      shift_file= false;
      logger << "No shift file";
    }
    else {
      int ret= fscanf(fp, "%e %e %e %e", &z, shift, shift+1, shift+2);
      assert(ret == 4);
      fclose(fp);
    }
  }
  mpi->bcast(&shift_file, 1, MPI_INT);  

  if(shift_file) {
    mpi->bcast(shift, 3, MPI_FLOAT);
    logger << "shift " 
	   << shift[0] << " " << shift[1] << " " << shift[2] << "\n";
  }

  if(read_xv_file) {
    read_pm_file_xv(filebase, redshift, mpi->index(), buffer_factor,
		    shift, particles);
    particles->boxsize = boxsize;
  }
  else {
    read_pm_file_zip(filebase, redshift, mpi->index(), buffer_factor, nc_node_dim,
      mesh_scale, shift,
      particles);
    boxsize= particles->boxsize;
  }
  //write_ascii(redshift, mpi->index(), particles);

  halos->boxsize= boxsize;

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
  const index_t nnode= tree->construct(particles, 
				 -buffer_factor*boxsize, 
				 (1.0f+buffer_factor)*boxsize);
  logger << "nnode " << nnode << "\n";


#ifdef FOFON
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
#endif

  //
  // Neighbor Search
  //
  
  mpi->print("neighbor search\n");
  logger.begin_timer(neighbors);
  KNeighbors knbrs;
  
  Particle* const p= particles->particle;
  const index_t np= particles->np_with_buffers;

#ifdef _OPENMP
#pragma omp parallel for private(knbrs)
#endif
  for(index_t i=0; i<np; ++i) {
    neighbor_search(p+i, &knbrs);
  }


#ifdef FOFON
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
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef HALO_EXCISE
  //
  // Spherical Overdensity Halo
  //
  logger.begin_timer(so_halo);
  /*sprintf(filename, "%s/%shalo%d.dat", 
    halo_dir, redshift, mpi->index());*/
  sprintf(filename, "%s%d/%shalo%d.dat",
	  filebase, mpi->index(), redshift, mpi->index());
  read_and_exchage_halo(filename, z, omega_m, mpi, halos, buffer_factor, shift);

  logger << "nhalo_local= " << halos->np_local << "\n";
  
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
    if(left < h[i].x[0] && h[i].x[0] < right &&
       left < h[i].x[1] && h[i].x[1] < right &&
       left < h[i].x[2] && h[i].x[2] < right) {
      int nph= mark_halo_particles(tree->root(), h[i].x, h[i].r);
      np_halo += nph;
      m_halo += h[i].mass/8.0f;
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
#else
  logger << "haloes not excised\n";
#endif


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
