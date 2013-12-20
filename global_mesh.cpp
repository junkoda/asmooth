#include <iostream>
#include <cmath>
#include "basic_types.h"
#include "mpi_interface.h"
#include "global_mesh.h"
#include "open.h"

using namespace std;

// TSC kernel
inline float f(const float x) {
  if(!(x >= 0.0f))
    cerr << "error x= " << x << endl; // debug!!!
  assert(x >= 0.0f);
  return x >= 1.0f ? 0.5f : x-0.5f*x*x;
}

inline float ker_integ(const float x1, const float x2)
{
  const float fac1= (float)(-1+2*(x1 > 0.0f));
  const float fac2= (float)(-1+2*(x2 > 0.0f));
  return fac2*f(fac2*x2) - fac1*f(fac1*x1);
}

void assign_on_global_mesh(Particles const * const particles,
			   MpiInterface const * const mpi,
			   float mesh[], const int nc, 
			   const int sign)
{
  const int nmesh= nc*nc*nc;
  float * const n= mesh;            // IGM particle number density
  float * const halo= n+nmesh;      // halo particle number density
  float * const c= halo + nmesh;    // IGM clumping
  float * const v= c + nmesh;

  Particle const * const p= particles->particle;
  const float boxsize= particles->boxsize;
  const float dx= boxsize/nc;
  const int nc2= nc*nc;
  
  const int np= particles->np_with_buffers;

  int const * const x_mpi= mpi->coordinate();
  float offset[]= {boxsize*x_mpi[0], boxsize*x_mpi[1], boxsize*x_mpi[2]};
  

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(index_t j=0; j<np; ++j) {
    float const * const y= p[j].x;
    if(y[0] < 0.0f || y[0] >= boxsize ||
       y[1] < 0.0f || y[1] >= boxsize ||
       y[2] < 0.0f || y[2] >= boxsize) {
      continue;
    }
       
    float x[3];
    x[0]= y[0] + offset[0];
    x[1]= y[1] + offset[1];
    x[2]= y[2] + offset[2];
    const float h= p[j].rk;
    assert(h > 0); // debug
    const float hinv= 1.0f/h;
    int left[3], right[3];
    for(int i=0; i<3; ++i) {
      left[i]=  (int)floor((x[i]-h)/boxsize * nc);
      right[i]= 1 + (int)ceil((x[i]+h)/boxsize * nc);
    }

    int ii[3];
    float s1;
    float ker[3];
    const float ds= dx*hinv;
    for(ii[2]= left[2]; ii[2]<right[2]; ii[2]++) {
      s1= (ii[2]*dx - x[2])*hinv;
      ker[2]= ker_integ(s1, s1+ds);
      for(ii[1]= left[1]; ii[1]<right[1]; ii[1]++) {
	s1= (ii[1]*dx - x[1])*hinv;
	ker[1]= ker_integ(s1, s1+ds);
	for(ii[0]= left[0]; ii[0]<right[0]; ii[0]++) {
	  s1= (ii[0]*dx - x[0])*hinv;
	  ker[0]= ker_integ(s1, s1+ds);
	  
	  int index= ((ii[2] + nc) % nc)*nc2 + 
	             ((ii[1] + nc) % nc)*nc + 
	             ((ii[0] + nc) % nc);
	  
	  const float ker3= ker[0]*ker[1]*ker[2];

	  if(sign*p[j].igrp <= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
	    n[index] += ker3;
#ifdef _OPENMP
#pragma omp atomic
#endif
	    c[index] += ker3*p[j].dens;
#ifdef _OPENMP
#pragma omp atomic
#endif
	    v[3*index] += ker3*p[j].v[0];
#ifdef _OPENMP
#pragma omp atomic
#endif
	    v[3*index+1] += ker3*p[j].v[1];
#ifdef _OPENMP
#pragma omp atomic
#endif
	    v[3*index+2] += ker3*p[j].v[2];
	  }
	  else if(sign*p[j].igrp > 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
	    halo[index] += ker3;
	  }  
	}
      }
    }
  }
}
