#include <iostream>
#include <cstdio>
#include "basic_types.h"
#include "file_reader.h"

using namespace std;

struct P3m_header
{
  int np_local;
  float a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, 
    cur_projection,cur_halofind,mass_p;
};

void read_pm_file_nbr(const char filename[],
		      const int dn[],
		      Particles* const particles,
		      const float buffer_factor);

//void read_pm_file(MpiInterface const * const mpi, const char filename[],
//		  Particles* const particles)
void read_pm_file(const char filename[], Particles* const particles)
{
  //
  // Read PM file and add particles
  //
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to read file: " << filename << endl;
    throw ErrorParticleReader();
  }
  
  //
  // Read header
  //
  P3m_header header;
  const int ret_read_header= fread(&header, sizeof(P3m_header), 1, fp);
  assert(ret_read_header == 1);

  if(particles->np_allocated < header.np_local) {
    cerr << "Not enough space for particles: "
         << particles->np_allocated << " < "
         << header.np_local << endl;
    throw ErrorParticleReader();
  }

  const index_t np= header.np_local;
  particles->np_local= np;
  particles->np_with_buffers= np;
  

  // get offset for global id
  //long np_local_l= (long) np;
  //long offset= 0;
  //MPI_Scan(&np_local_l, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  //particles->id_begin= offset - np + 1;  // id starts at 1



  //
  // Read particles
  //
  //particle_id_t id= particles->id_begin;
  Particle* p= particles->particle;
  float xv[6];
  for(int j=0; j<np; ++j) {
    const int ret= fread(xv, sizeof(float), 6, fp);
    assert(ret == 6);
    p->x[0]= xv[0]; p->x[1]= xv[1]; p->x[2]= xv[2];
    p->v[0]= xv[3]; p->v[1]= xv[4]; p->v[2]= xv[5];
    p->rk= 0.0f;
    p->dens= 0.0f;
    //p->id= id++;
    p->igrp= 0;

    ++p;
  }

  //particles->particle[np].id = -1;

  const int ret_fclose= fclose(fp);
  assert(ret_fclose == 0);

}

void read_pm_file_nbr(const char filename[],
		      const int dn[],
		      Particles* const particles,
		      const float buffer_factor,
		      const float shift[]);


void read_pm_file2(const char filename[], Particles* const particles,
		   const float shift[])
{
  //
  // Read PM file and add particles (all at once)
  //
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to read file: " << filename << endl;
    throw ErrorParticleReader();
  }
  
  //
  // Read header
  //
  P3m_header header;
  const int ret_read_header= fread(&header, sizeof(P3m_header), 1, fp);
  assert(ret_read_header == 1);

  if(particles->np_allocated < header.np_local) {
    cerr << "Not enough space for particles: "
         << particles->np_allocated << " < "
         << header.np_local << endl;
    throw ErrorParticleReader();
  }

  const index_t np= header.np_local;
  particles->np_local= np;
  particles->np_with_buffers= np;
  

  float* const buf= (float*) particles->particle;
  float* xv= buf + (sizeof(Particle)/sizeof(float)*particles->np_allocated -
                   6*(np+1));

  // get offset for global id
  //long np_local_l= (long) np;
  //long offset= 0;
  //MPI_Scan(&np_local_l, &offset, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  //particles->id_begin= offset - np + 1;  // id starts at 1



  //
  // Read particles
  //
  //particle_id_t id= particles->id_begin;
  const int ret= fread(xv, sizeof(float), 6*np, fp);
  assert(ret == 6*np);

  const int ret_fclose= fclose(fp);
  assert(ret_fclose == 0);

  //Put xv data into Particle format
  Particle* p= particles->particle;

  for(int j=0; j<np; ++j) {
    p->x[0]= xv[0] - shift[0];
    p->x[1]= xv[1] - shift[1];
    p->x[2]= xv[2] - shift[2];
    p->v[0]= xv[3]; p->v[1]= xv[4]; p->v[2]= xv[5];
    p->rk= 0.0f;
    p->dens= 0.0f;
    //p->id= id++;
    p->igrp= 0;

    ++p;
    xv += 6;
  }

  //particles->particle[np].id = -1;
}


//
// Read particles from 27 nbr files
//
void read_pm_file_all(const char redshift[], 
		      const int file_index, const int file_nc,
		      Particles* const particles,
		      const float buffer_factor,
		      const float shift[])
{
  char filename[128];

  sprintf(filename, "%sxv%d.dat", redshift, file_index);
  read_pm_file(filename, particles);

  int n[3];
  n[2]= file_index / (file_nc*file_nc); 
  n[1]= (file_index - file_nc*file_nc*n[2])/file_nc; 
  n[0]= file_index % file_nc;
  assert(file_nc*file_nc*n[2] + file_nc*n[1] + n[0] == file_index);
  
  // read neighbor file
  int dn[3];
  for(dn[2]= -1; dn[2]<2; ++dn[2]) {
   for(dn[1]= -1; dn[1]<2; ++dn[1]) {
    for(dn[0]= -1; dn[0]<2; ++dn[0]) {
      int index=   ((n[2]+dn[2]+file_nc) % file_nc)*file_nc*file_nc
	         + ((n[1]+dn[1]+file_nc) % file_nc)*file_nc
	         + ((n[0]+dn[0]+file_nc) % file_nc);
      //cerr << index << endl; // debug
      sprintf(filename, "%sxv%d.dat", redshift, index);
      if(index != file_index)       
	read_pm_file_nbr(filename, dn, particles, buffer_factor, shift);
    }
   }
  }
}

void read_pm_file_nbr(const char filename[],
		      const int dn[],
		      Particles* const particles,
		      const float buffer_factor,
		      const float shift[])
{
  //
  // Read PM file and add neighbor particles
  //
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to read file: " << filename << endl;
    throw ErrorParticleReader();
  }
  
  //
  // Read header
  //
  P3m_header header;
  const int ret_read_header= fread(&header, sizeof(P3m_header), 1, fp);
  assert(ret_read_header == 1);

  const index_t np= header.np_local;
  const float boxsize= particles->boxsize;
  const float dx_buffer= buffer_factor*boxsize;


  //
  // Read particles
  //
  Particle* p= particles->particle + particles->np_with_buffers;
  Particle* const p_end= particles->particle + particles->np_allocated;
  index_t np_with_buffers= particles->np_with_buffers;

  float xv[6];
  for(int j=0; j<np; ++j) {
    const int ret= fread(xv, sizeof(float), 6, fp);
    assert(ret == 6);
    if(dn[0]*xv[0] < -((1-dn[0])/2)*boxsize + dx_buffer &&
       dn[1]*xv[1] < -((1-dn[1])/2)*boxsize + dx_buffer &&
       dn[2]*xv[2] < -((1-dn[2])/2)*boxsize + dx_buffer) {
     
      //assert(p != p_end);
      if(p == p_end) {
	cerr << "No enough space for buffers: "
             << filename << " " 
	     << particles->np_local << " "
	     << particles->np_allocated << endl;
      }

      p->x[0]= xv[0] + dn[0]*boxsize - shift[0];
      p->x[1]= xv[1] + dn[1]*boxsize - shift[1];
      p->x[2]= xv[2] + dn[2]*boxsize - shift[2];
      p->v[0]= xv[3]; p->v[1]= xv[4]; p->v[2]= xv[5];
      p->rk= 0.0f;
      p->dens= 0.0f;
      p->igrp= 0;
      
      ++p;
      ++np_with_buffers;
    }
  }

  particles->np_with_buffers= np_with_buffers;

  const int ret_fclose= fclose(fp);
  assert(ret_fclose == 0);
}

