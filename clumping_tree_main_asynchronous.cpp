//
// Calculate density, clumping, velocity field from adaptive kernel
//   asynchronous version: each node reads all 27 nbr files without communication
//


#include <iostream>
#include <deque>
#include <mpi.h>
#include "option.h"
#include "logger.h"
#include "basic_types.h"
// #include "mpi_interface.h"
#include "file_reader.h"
#include "exchange_buffer.h"
#include "kdtree_balanced.h"
#include "nbr_search.h"
#include "fof.h"
#include "coarse_mesh.h"
#include "halo_file.h"

using namespace std;

void for_each_file(const int file_index,
		   const int file_nc,
		   const char redshift[],
		   const char particle_dir[],
		   const char halo_dir[],
		   Logger& logger,
		   Particles* const particles,
		   ParticleSet<Halo>* halos,
		   const float buffer_factor,
		   KDTree* const tree,
		   int* const temp_buffer,
		   float* const mesh,
		   deque<int>& ncs);

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
  //
  // MPI setup
  //
  MPI_Init(&argc, &argv);
  int my_rank, n_node;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_node);
  if(n_node == 1) {
    cerr << "Need more than 1 node\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  //MpiInterface* const mpi= new MpiInterface();


  if(argc == 1 && my_rank == 0) {
    cout << "clumping_tree [options]\n"
	 << "\t-pm_redshift <z>; p3m file <z>xv<i>.dat\n"
	 << "\t-pdir <particle directory>; -pdir particle_data\n"
	 << "\t-allocate <MB for particles>\n"
	 << "\t-boxsize <boxsize>\n"
	 << "\t-buffer_factor 0.05\n"
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
  op.set_default("-nc", "256");       // nc^3 is number of mesh in this node
  op.set_default("-nc_file", "8"); // nc_file^3 = number of files per redshift
  //op.set_default("-ncf", "0"); // finemesh clumping

  const float boxsize= op.get_float("-boxsize");
  assert(boxsize > 0.0f);
  const float buffer_factor= op.get_float("-buffer_factor");
  assert(buffer_factor > 0.0f);
  //const int nc= op.get_int("-nc");
  const int nc_file= op.get_int("-nc_file");
  const int n_file= nc_file*nc_file*nc_file;

  logger << "boxsize " << boxsize << "\n";
  logger << "buffer_factor " << buffer_factor << "\n";
  
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
  
  

  //
  // open redshift list
  //
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
    MPI_Abort(MPI_COMM_WORLD, 1);
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
  // For each snapshot
  //
  //for(deque<string>::iterator z= redshifts.begin();
  //z != redshifts.end(); ++z) {
  //bool finished= false;
  //while(finished) {
  // not done

  const int nmax= n_file*redshifts.size(); //debug! nc_file^3*n_redshift
  if(my_rank == 0) {
    //
    // 0: distribute integer (only)
    //
    int count= 0;
    int recv= 0;
    int done= 1;
    MPI_Status status;

    while(1) {
      MPI_Recv(&recv, 1, MPI_INT, 
	       MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if(recv == 1) {
	done++;
	if(done == n_node) break; // all nodes are done
      }
      MPI_Send(&count, 1, MPI_INT,
	       status.MPI_SOURCE, 1, MPI_COMM_WORLD);
      count++;
    }
  }
  else {
    //
    // id > 0: get number and process
    //
    int dummy= 0;
    int get= 0;
    MPI_Status status;
    while(1) {
      MPI_Send(&dummy, 1, MPI_INT,
	       0, 1, MPI_COMM_WORLD);
      MPI_Recv(&get, 1, MPI_INT,
	       0, 1, MPI_COMM_WORLD, &status);
      //cout << get << endl;
      if(get < nmax) {
	printf("file %3d by %d\n", get, my_rank);

	//work
	const int n_redshift= get / n_file;
	const int file_index= get % n_file;

	for_each_file(file_index,
		      nc_file,
		      redshifts[n_redshift].c_str(), // redshift
		      op["-pdir"].c_str(), // particle dir
		      op["-hdir"].c_str(), // halo dir
		      logger,
		      particles,
		      halos,
		      buffer_factor,
		      tree,
		      temp,
		      mesh,
		      ncs);
      }
      else {
	int recv_done= 1;
	MPI_Send(&recv_done, 1, MPI_INT,
		 0, 1, MPI_COMM_WORLD);
	break;
      }
    }
  }

  logger << "done all\n";

  delete tree;
  //delete mpi;
  free(mesh);
  free(temp);

  MPI_Finalize();

  return 0;
}

void for_each_file(const int file_index,
		   const int file_nc,
		   const char redshift[],
		   const char particle_dir[],
		   const char halo_dir[],
		   Logger& logger,
		   Particles* const particles,
		   ParticleSet<Halo>* halos,
		   const float buffer_factor,
		   KDTree* const tree,
		   int* const temp_buffer,
		   float* const mesh,
		   deque<int>& ncs)
{
  const float boxsize= particles->boxsize;
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  //
  // Read File
  //
  if(my_rank == 1)
    printf("reading particles\n");
  logger.begin_timer(io);

  logger << "z " << redshift << "\n";
  char filename[256];

  float shift[]= {0.0f, 0.0f, 0.0f};
  // read all particles from file
  sprintf(filename, "%s/%s", particle_dir, redshift);
  read_pm_file_all(filename, 
		   file_index, file_nc,
		   particles,
		   buffer_factor,
		   shift);

  //sprintf(filename, "%s/%sxv%d.dat", 
  //particle_dir, redshift, mpi->index());
  //read_pm_file(filename, particles);

  if(particles->np_local <= 0) {
    cerr << "No particles\n";
    MPI_Abort(MPI_COMM_WORLD, 1);//mpi->abort("No particles\n");
  }

  logger << "np_local " << particles->np_local << "\n";
  logger << "np_with_buffers " << particles->np_with_buffers << "\n";

  //
  // Construct KDTree
  //
  if(my_rank == 1)
    printf("constructing kdtree\n");
    //mpi->print("constructing kdtree\n");
  logger.begin_timer(construction);
  //KDTree* tree= new KDTreeBalanced();
  const index_t nnode= tree->construct(particles, 
				 -buffer_factor*boxsize, 
				 (1.0f+buffer_factor)*boxsize);
  logger << "nnode " << nnode << "\n";


  //
  // FOF
  //
  if(my_rank == 1)
    printf("fof\n");
  //mpi->print("fof\n");
  logger.begin_timer(fof);
  const index_t ngrp= group_fof(particles, tree, 0.4);
  logger << "ngrp " << ngrp << "\n";
  const index_t nhalo= remove_small_groups(particles, temp_buffer, ngrp);
  logger << "nhalo " << nhalo << "\n";


  //
  // Neighbor Search
  //
  
  if(my_rank == 1)
    printf("neighbor search\n");
  //mpi->print("neighbor search\n");
  logger.begin_timer(neighbors);
  KNeighbors knbrs;
  
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
  logger.begin_timer(mesh_assign);

  for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
    const int nc= *p; assert(nc > 0);
    clear_mesh(mesh, nc);
    assign_on_mesh(particles, mesh, nc); // local

    //logger.begin_timer(output);
    sprintf(filename, "fof/nc%d/%s", nc, redshift);
    write_mesh_separate(filename, file_index, mesh, nc, boxsize);
  }

  //logger.begin_timer(imbalance);

  //
  // Spherical Overdensity Halo
  //
  logger.begin_timer(so_halo);
  /*
  sprintf(filename, "%s/%shalo%d.dat", 
	  halo_dir, redshift, mpi->index());
  read_and_exchage_halo(filename, mpi, halos, buffer_factor, full_box);
  */
  sprintf(filename, "%s/%s", halo_dir, redshift);
  read_halo_file_all(filename, 
		     file_index, file_nc,
		     halos, buffer_factor, shift);

  const Halo* h= halos->particle;
  const index_t n_sohalo= halos->np_with_buffers;

  // Mark Halo Particles
  int np_halo= 0;
  double m_halo= 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:np_halo, m_halo)
#endif
  for(index_t i=0; i<n_sohalo; ++i) {
    np_halo += mark_halo_particles(tree->root(), h[i].x, h[i].r);
    m_halo += h[i].mass/8.0;
  }

  logger << "so_halo_total " << np_halo << " " << m_halo << "\n";

  //  
  // Mesh Assignment (SO)
  //
  logger.begin_timer(mesh_assign);

  for(deque<int>::iterator p= ncs.begin(); p != ncs.end(); ++p) {
    const int nc= *p;
    clear_mesh(mesh, nc);

    assign_on_mesh(particles, mesh, nc, -1);

    sprintf(filename, "so/nc%d/%s", nc, redshift);
    write_mesh_separate(filename, file_index, mesh, nc, boxsize);
  }
  //logger.begin_timer(output);
  
  //
  // Done
  //

  //logger.begin_timer(imbalance);

  logger.end_timer();
    
  if(my_rank == 1) {
    logger.print_log(cout);
    printf("bye\n");
  }
  //mpi->print("bye\n");
}
