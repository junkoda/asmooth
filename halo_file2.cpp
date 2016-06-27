
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "kdtree_balanced.h"
#include "halo_file.h"
#include "exchange_buffer_template.h"
#include "nbr_search.h"
#include <stdint.h>
using namespace std;

static const int halo_ndat= 22; // number of float per halo (old ver. 17)
static const int halo_x= 0;     // offset for halo position (starting 0)
static const int halo_m= 4;     // offset for halo mass
static const int halo_r= 6;     // offset for halo radius
static const float halo_delta= 200.0f; // overdensity

void correct_halo_coordinate(ParticleSet<Halo>* const halos);
void correct_halo_coordinate_full(ParticleSet<Halo>* const halos,
				  const int icpu[]);
int read_halo_file_nbr(const char filename[], 
		       const int nc_node[],
		       ParticleSet<Halo>* const halos,
		       const float buffer_factor,
		       const float shift[]);

// Read halo_data/*halo<n>.dat
// 0-2  : hpos(3)
// 3    : mass_vir
// 4    : mass_odc
// 5    : r_vir
// 6    : r_odc
// ...



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


  //const float delta= 745.6046565f; // 4pi/3*178
  const float boxsize= halos->boxsize;
  int n= 0;
  float buf[halo_ndat];
  Halo* halo= halos->particle;

  int header[10];
  fread(header, sizeof(int), 3, fp); // header??

  while(fread(buf, sizeof(float), halo_ndat, fp) == halo_ndat) {
    halo->x[0]= fmod(buf[halo_x], boxsize) - shift[0];
    halo->x[1]= fmod(buf[halo_x + 1], boxsize) - shift[1];
    halo->x[2]= fmod(buf[halo_x + 2], boxsize) - shift[2];
    const float r= halo->r= buf[halo_r];
    const float m= buf[halo_m];
    halo->mass= m;
    assert(fabs(m/(4.0f/3.0f*M_PI*r*r*r) - halo_delta) < 10.0f);

#ifdef PID_FLAG
    // 50 8-byte integer for id and 6 4-byte floats for xv
    fseek(fp, 50*(8+4*6), SEEK_CUR);
#endif
    
    ++n;
    ++halo;
  }

  halos->np_local= n;
  halos->np_with_buffers= n;

  return fclose(fp) == 0;
}

