
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include <stdint.h>
#include "kdtree_balanced.h"
#include "halo_file.h"
#include "exchange_buffer_template.h"
#include "nbr_search.h"
#include "endian.h"

using namespace std;

// Halo File format
#ifdef M200
static const int halo_m= 4;     // offset for halo mass
static const int halo_r= 6;     // offset for halo radius
static const float delta200= 200.0f; // overdensity
static const float torr= 10.0f;
#else
static const int halo_m= 3;     // offset for halo mass (Bryan and Norman)
static const int halo_r= 5;     // offset for halo radius
static const float torr= 10.0f;
#endif

static const int halo_ndat= 22; // number of float per halo (old ver. 17)
static const int halo_x= 0;     // offset for halo position (starting 0)
static const int npid= 50;

void correct_halo_coordinate(ParticleSet<Halo>* const halos);
void correct_halo_coordinate_full(ParticleSet<Halo>* const halos,
				  const int icpu[]);
int read_halo_file(const char filename[], ParticleSet<Halo>* const halos,
		   const float shift[], const float halo_delta);

// Read halo_data/*halo<n>.dat
// 0-2  : hpos(3)
// 3    : mass_vir
// 4    : mass_odc
// 5    : r_vir
// 6    : r_odc
// ...


double halo_delta_vir(const double z, const double omega_m)
{
#ifdef M200
  return delta200;
#else
  assert(omega_m > 0.0);
  const double a= 1.0/(1.0 + z);
  const double omega_l= 1.0 - omega_m;
  const double omega_z= omega_m/(omega_m + omega_l*a*a*a);
  const double x= omega_z - 1.0;
  
  return (18.0*M_PI*M_PI + 82.0*x - 39.0*x*x)*
    ((omega_m + omega_l*a*a*a)/(omega_m));
#endif
}

void read_and_exchage_halo(const char filename[],
			   const float z, const float omega_m,
			   MpiInterface const * const mpi,
			   ParticleSet<Halo>* const halos,
			   const float buffer_factor,
			   const float shift[]
			   )
{
  float delta_halo = halo_delta_vir(z, omega_m);
  if(mpi->rank() == 0) {
    printf("omega_m = %.4f\n", omega_m);
    printf("halo_delta = %.3f at z= %.4f\n", delta_halo, z);
  }

  read_halo_file(filename, halos, shift, delta_halo);
  /*
  if(full_box)
    correct_halo_coordinate_full(halos, mpi->coordinate());
  else
    correct_halo_coordinate(halos);
  */

  exchange_buffer_template(mpi, halos, buffer_factor);
}

int read_halo_file(const char filename[], ParticleSet<Halo>* const halos,
		const float shift[], const float halo_delta)
{
  int nhalo=0;
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open halo data: " << filename << endl;
    assert(false);
  }


  const float boxsize= halos->boxsize;
  int n= 0;
  float buf[halo_ndat];
  Halo* halo= halos->particle;

  int header[10];
  fread(header, sizeof(int), 3, fp); // header??

  while(fread(buf, sizeof(float), halo_ndat, fp) == halo_ndat) {
    halo->x[0]= fmod(endian_float(buf[halo_x]), boxsize) - shift[0];
    halo->x[1]= fmod(endian_float(buf[halo_x + 1]), boxsize) - shift[1];
    halo->x[2]= fmod(endian_float(buf[halo_x + 2]), boxsize) - shift[2];
    const float r= halo->r= endian_float(buf[halo_r]);
    const float m= endian_float(buf[halo_m]);
    halo->mass= m;
    //assert(fabs(m/(4.0f/3.0f*M_PI*r*r*r) - halo_delta) < torr);
    if(fabs(m/(4.0f/3.0f*M_PI*r*r*r) - halo_delta) >= torr) {
      cerr << "Halo over density in file seems incorrect\n";
      cerr << "halo_delta " << halo_delta << endl;
      cerr << "delta data = " << m/(4.0f/3.0f*M_PI*r*r*r) << endl;
      cerr << "Do not agree with in torr = " << torr << endl;
      assert(false);
    }

#ifdef PID_FLAG
    // 50 8-byte integer for id and 6 4-byte floats for xv
    int seek_ret= fseek(fp, npid*(8+4*6), SEEK_CUR);
    assert(seek_ret == 0);
#endif
    
    ++n;
    ++halo;
  }

  halos->np_local= n;
  halos->np_with_buffers= n;

  int fret= fclose(fp); assert(fret == 0);
  return fret == 0;
}

