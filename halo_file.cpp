
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "kdtree_balanced.h"
#include "halo_file.h"
#include "exchange_buffer_template.h"
#include "nbr_search.h"

using namespace std;

static const int ndat= 17; // number of float per halo (old ver. 17)

void correct_halo_coordinate(ParticleSet<Halo>* const halos);
void correct_halo_coordinate_full(ParticleSet<Halo>* const halos,
				  const int icpu[]);
int read_halo_file_nbr(const char filename[], 
		       const int nc_node[],
		       ParticleSet<Halo>* const halos,
		       const float buffer_factor,
		       const float shift[]);

// Read halo_data/*halo<n>.dat
// See makesourcelist.f90 for data format
// Column 1-3:   halo_pos(:)
// Column 4-6:   x_mean(:)
// Column 7-9:   v_mean(:)
// Column 10-12: l(:)
// Column 13 :   v_disp
// Column 14 :   radius_calc
// Column 15 :   halo_mass
// Column 16 :   imass
// Column 17 :   halo_mass_no_corr


void read_and_exchage_halo(const char filename[],
			   MpiInterface const * const mpi,
			   ParticleSet<Halo>* const halos,
			   const float buffer_factor,
			   const float shift[])
{
  read_halo_file(filename, halos, shift);

  /*
  if(full_box)
    correct_halo_coordinate_full(halos, mpi->coordinate());
  else
    correct_halo_coordinate(halos);
  */

  exchange_buffer_template(mpi, halos, buffer_factor);
}

int read_halo_file(const char filename[], ParticleSet<Halo>* const halos,
		   const float shift[])
{
  int nhalo=0;
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open halo data: " << filename << endl;
    assert(false);
  }
  
  fread(&nhalo, sizeof(int), 1, fp);

  //const float delta= 745.6046565f; // 4pi/3*178
  const float boxsize= halos->boxsize;
  int n= 0;
  float buf[ndat];
  Halo* halo= halos->particle;

  while(fread(buf, sizeof(float), ndat, fp) != 0) {
    halo->x[0]= fmod(buf[0], boxsize) - shift[0];
    halo->x[1]= fmod(buf[1], boxsize) - shift[1];
    halo->x[2]= fmod(buf[2], boxsize) - shift[2];
    const float r= halo->r= buf[13];
    //const float m= halo->mass= buf[14];
    const float m= buf[14];
    halo->mass= buf[15];
    //halo->r= pow(buf[14]/delta, 1.0f/3); // based on gridmass (correct)
    //halo->r= buf[13]; // (imass/delta)^(1/3)
    printf("%e %e %e %e %e\n", halo->x[0], halo->x[1], halo->x[2], m, r);
    //assert(abs(m/(4.0f/3.0f*M_PI*r*r*r) - 178.0f) < 0.1);

    ++n;
    ++halo;
  }

  halos->np_local= n;
  halos->np_with_buffers= n;
  //cerr << "nhalo_so_local " << n << endl;

  return fclose(fp) == 0;
}

/*
void remove_halo_particles(vector<Halo> const * const pv,
			   TreeSystem_for_knbr<ParticleData>* const tree)
{
  for(vector<Halo>::const_iterator p= pv->begin(); p != pv->end(); ++p)
    tree->remove_r_nbr(p->x, p->r);
}
*/

// not using anymore (corrected while reading)
/*
void correct_halo_coordinate(ParticleSet<Halo>* const halos)
{
  Halo * const h= halos->particle;
  const float boxsize= halos->boxsize;
  const index_t np= halos->np_local;
  
  for(index_t i=0; i<np; ++i) {
    h[i].x[0]= fmod(h[i].x[0], boxsize);
    h[i].x[1]= fmod(h[i].x[1], boxsize);
    h[i].x[2]= fmod(h[i].x[2], boxsize);
  }
}

void correct_halo_coordinate_full(ParticleSet<Halo>* const halos,
				  const int icpu[])
{
  Halo * const h= halos->particle;
  const float boxsize= halos->boxsize;
  const index_t nhalo= halos->np_local;
  
  //cerr << "icpu " << icpu[0] << " " << icpu[1] << " " << icpu[2] << endl;

  const float dx= 0.05f*boxsize;

  for(int i=0; i<nhalo; ++i) {
    h[i].x[0] -= icpu[0]*boxsize;
    h[i].x[1] -= icpu[1]*boxsize;
    h[i].x[2] -= icpu[2]*boxsize;
    
    assert(-dx < h[i].x[0] && h[i].x[0] < boxsize+dx &&
	   -dx < h[i].x[1] && h[i].x[1] < boxsize+dx &&
	   -dx < h[i].x[2] && h[i].x[2] < boxsize+dx);
  }
}
*/

//
// Read all neighbor halos
//
int read_halo_file_all(const char redshift[], 
		       const int file_index, const int file_nc,
		       ParticleSet<Halo>* const halos,
		       const float buffer_factor,
		       const float shift[])
{
  int n_file[3];
  n_file[2]= file_index / (file_nc*file_nc); 
  n_file[1]= (file_index - file_nc*file_nc*n_file[2])/file_nc; 
  n_file[0]= file_index % file_nc;
  assert(file_nc*file_nc*n_file[2] + file_nc*n_file[1] + n_file[0] 
	 == file_index);

  char filename[128];
  sprintf(filename, "%shalo%d.dat", redshift, file_index);
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open halo data: " << filename << endl;
    assert(false);
  }

  int nhalo=0;
  fread(&nhalo, sizeof(int), 1, fp);

  const float boxsize= halos->boxsize;
  //const float delta= 745.6046565f; // 4pi/3*178
  int n= 0;
  float buf[ndat];
  Halo* halo= halos->particle;

  while(fread(buf, sizeof(float), ndat, fp) != 0) {
    halo->x[0]= buf[0] - n_file[0]*boxsize - shift[0];
    halo->x[1]= buf[1] - n_file[1]*boxsize - shift[1];
    halo->x[2]= buf[2] - n_file[2]*boxsize - shift[2];
    //halo->mass= buf[14];
    //halo->r= pow(buf[14]/delta, 1.0f/3); // based on gridmass (correct)
    //halo->r= buf[13]; // this should be OK, too
    const float r= halo->r= buf[13];
    const float m= halo->mass= buf[14];
    assert(abs(m/(4.0f/3.0f*M_PI*r*r*r) - 178.0f) < 0.1);

    ++n;
    ++halo;
  }

  halos->np_local= n;
  halos->np_with_buffers= n;

  //
  // read neighbor halos
  //
  int dn[3];
  for(dn[2]= -1; dn[2]<2; ++dn[2]) {
   for(dn[1]= -1; dn[1]<2; ++dn[1]) {
    for(dn[0]= -1; dn[0]<2; ++dn[0]) {
      int index=   ((n_file[2]+dn[2]+file_nc) % file_nc)*file_nc*file_nc
	         + ((n_file[1]+dn[1]+file_nc) % file_nc)*file_nc
	         + ((n_file[0]+dn[0]+file_nc) % file_nc);
      sprintf(filename, "%shalo%d.dat", redshift, index);
      if(index != file_index)       
	read_halo_file_nbr(filename, n_file, halos, buffer_factor, shift);
    }
   }
  }
  return 0;
}

int read_halo_file_nbr(const char filename[], 
		       const int n_file[],
		       ParticleSet<Halo>* const halos,
		       const float buffer_factor,
		       const float shift[])
{
  int nhalo=0;
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open halo data: " << filename << endl;
    assert(false);
  }
  
  fread(&nhalo, sizeof(int), 1, fp);

  const float boxsize= halos->boxsize;
  const float left= -buffer_factor*boxsize;
  const float right= (1.0f+buffer_factor)*boxsize;

  const float delta= 745.6046565f; // 4pi/3*178
  int n= 0;
  float buf[ndat];
  Halo* halo= halos->particle + halos->np_with_buffers;
  Halo* const halo_end= halos->particle + halos->np_allocated;

  while(fread(buf, sizeof(float), ndat, fp) != 0) {
    halo->x[0]= buf[0] - n_file[0]*boxsize - shift[0];
    halo->x[1]= buf[1] - n_file[1]*boxsize - shift[1];
    halo->x[2]= buf[2] - n_file[2]*boxsize - shift[2];
    //halo->mass= buf[14];
    //halo->r= pow(buf[14]/delta, 1.0f/3); // based on gridmass
    //halo->r= buf[13]; // (imass/delta)^(1/3) // probably same

    const float r= halo->r= buf[13];
    const float m= halo->mass= buf[14];
    assert(abs(m/(4.0f/3.0f*M_PI*r*r*r) - 178.0f) < 0.1);

    if(left < halo->x[0] && halo->x[0] < right &&
       left < halo->x[1] && halo->x[1] < right &&
       left < halo->x[2] && halo->x[2] < right) {
      assert(halo != halo_end);
      ++n;
      ++halo;
    }
  }

  halos->np_with_buffers += n;

  return fclose(fp) == 0;
}

