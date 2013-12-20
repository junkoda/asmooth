#include <iostream>
#include <cmath>
#include "basic_types.h"
#include "coarse_mesh.h"
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

void clear_mesh(float mesh[], const int nc)
{
  const int nmesh= nc*nc*nc*6;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j=0; j<nmesh; ++j)
    mesh[j]= 0.0f;
}

double mesh_total(const float mesh[], const int nc)
{
  const int nmesh= nc*nc*nc;
  double mtot= 0.0;
  for(int j=0; j<nmesh; ++j)
    mtot += mesh[j];

  return mtot;
}

void assign_on_mesh(Particles const * const particles,
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

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(index_t j=0; j<np; ++j) {
    float const * const x= p[j].x;
    const float h= p[j].rk;
    assert(h > 0); // debug
    const float hinv= 1.0f/h;
    //const int offset= nmesh*(sign*p[j].igrp > 0); 
    int left[3], right[3];
    for(int i=0; i<3; ++i) {
      left[i]=  max(0,  (int)((x[i]-h)/boxsize * nc));
      right[i]= min(nc, 1+(int)((x[i]+h)/boxsize * nc));
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
	  
	  int index= ii[2]*nc2 + ii[1]*nc + ii[0];
	  
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

void write_mesh(char filename[], float* mesh, const int nc, const float boxsize)
{
  //
  // ** output **
  // int nc
  // int dof (number of floats per mesh)
  // [n, nhalo, sum(n), n*v[3]]
  FILE* fp= open(filename, "w");
  assert(fp);
  const int dof= 6;
  fwrite(&nc, sizeof(int), 1, fp);  // number of mesh per dim
  fwrite(&dof, sizeof(int), 1, fp); // number of float per mesh
  fwrite(&boxsize, sizeof(float), 1, fp);
  fwrite(mesh, sizeof(float), dof*nc*nc*nc, fp);
                                    // n halo n*c n*v[3]
  fclose(fp);
}

void write_mesh_separate(const char base_filename[], const int inode, 
			 float *const mesh, const int nc, const float boxsize)
{
  const int nmesh= nc*nc*nc;
  const float vol= pow(boxsize/nc, 3.0f);

  float * const n= mesh;
  float * const halo= mesh+nmesh;
  float * const c= halo+nmesh;
  float * const v= c+nmesh;

  // normalize c, v, n
  for(int j=0; j<nmesh; ++j) {
    if(n[j] > 0.0f) {
      c[j] *= vol/(n[j]*n[j]);
      v[3*j]   /= n[j];
      v[3*j+1] /= n[j];
      v[3*j+2] /= n[j];
    }
    else {
      c[j]= v[3*j]= v[3*j+1]= v[3*j+2]= 0.0f;
    }
    n[j] *= 8.0f;      // normalization factor (particle mass= 8)
    halo[j] *= 8.0f;
  }
  
  char filename[128];
  
  // write 'n'
  sprintf(filename, "%sn%d.dat", base_filename, inode);
  FILE* fp= open(filename, "w"); assert(fp);
  fwrite(&nc, sizeof(int), 1, fp);
  fwrite(n, sizeof(float), nmesh, fp);
  fclose(fp);

  // write 'halo'
  sprintf(filename, "%sh%d.dat", base_filename, inode);
  fp= open(filename, "w"); assert(fp);
  fwrite(&nc, sizeof(int), 1, fp);
  fwrite(halo, sizeof(float), nmesh, fp);
  fclose(fp);

  // write 'c'
  sprintf(filename, "%sc%d.dat", base_filename, inode);
  fp= open(filename, "w"); assert(fp);
  fwrite(&nc, sizeof(int), 1, fp);
  fwrite(c, sizeof(float), nmesh, fp);
  fclose(fp);

  // write 'vel'
  sprintf(filename, "%svel%d.dat", base_filename, inode);
  fp= open(filename, "w"); assert(fp);
  fwrite(&nc, sizeof(int), 1, fp);
  fwrite(v, sizeof(float), 3*nmesh, fp);
  fclose(fp);


  // Units
  // n ~ 8*sum(1)
  // c ~ [sum(n_i)/V] / [sum(1)/V]^2   # dimmensionless
  // vel ~ sum(v)/sum(1)               # mean velocity
}



void assign_on_density_mesh(Particles const * const particles,
			 float mesh[], const int nc,
			 const int sign)
{
  //const int nmesh= nc*nc*nc;
  float * const n= mesh;            // IGM particle number density

  Particle const * const p= particles->particle;
  const float boxsize= particles->boxsize;
  const float dx= boxsize/nc;
  const int nc2= nc*nc;
  
  const int np= particles->np_with_buffers;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(index_t j=0; j<np; ++j) {
    float const * const x= p[j].x;
    const float h= p[j].rk;
    assert(h > 0); // debug
    const float hinv= 1.0f/h;
    //const int offset= nmesh*(sign*p[j].igrp > 0); 
    int left[3], right[3];
    for(int i=0; i<3; ++i) {
      left[i]=  max(0,  (int)((x[i]-h)/boxsize * nc));
      right[i]= min(nc, 1+(int)((x[i]+h)/boxsize * nc));
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
	  
	  int index= ii[2]*nc2 + ii[1]*nc + ii[0];
	  
	  const float ker3= ker[0]*ker[1]*ker[2];

	  if(sign*p[j].igrp <= 0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
	    n[index] += ker3;
	  }
	}
      }
    }
  }
}
