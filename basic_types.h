#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H 1

#include <cstdlib>
#include <assert.h>

typedef int index_t; // Must be large enough to label np_local (np this node)
typedef long particle_id_t;   // Must be large enough to lable np_total

struct Particle {
  float x[3], v[3];
  float dens, rk;
  //particle_id_t id;
  int igrp;
};

struct Halo {
  float x[3];
  float r, mass;
};


//
// Particles: A set of particles
//
struct Particles {
  Particles();
  ~Particles();

  Particle* particle;
  index_t np_allocated, np_local, np_with_buffers;
  particle_id_t id_begin, np_total;
  float boxsize;

  void allocate(const index_t n_alloc);
};

/*
struct Halos {
  Halos();
  ~Halos();
  Halo* halo;
  int np_allocated, np_local, np_with_buffers;
  float boxsize;

  void allocate(const int n_alloc);
}
*/

template<class T>
struct ParticleSet {
  ParticleSet();
  ~ParticleSet();
  T* particle;
  index_t np_allocated, np_local, np_with_buffers;
  float boxsize;

  void allocate(const index_t n_alloc);
  void allocate(void* buf, const index_t n_alloc);
};

template<class T>
ParticleSet<T>::ParticleSet() : 
  particle(0), np_allocated(0)
{

}

template<class T>
ParticleSet<T>::~ParticleSet()
{
  //free(particle); particle= 0;
}

template<class T>
void ParticleSet<T>::allocate(const int n_alloc)
{
  free(particle);
  particle= (T*) malloc(sizeof(T)*n_alloc);
  assert(particle);
  np_allocated= n_alloc;
  np_local= np_with_buffers= 0;
}

template<class T>
void ParticleSet<T>::allocate(void *buf, const index_t n_alloc)
{
  particle= (T*) buf;
  assert(particle);
  np_allocated= n_alloc;
  np_local= np_with_buffers= 0;
}


#endif
