#include "basic_types.h"
#include <cstdlib>
#include <assert.h>

//
// Particles: A Set of Particles
//
Particles::Particles() :
  particle(0), 
  np_allocated(0), np_local(0), np_with_buffers(0),
  id_begin(0), np_total(0),
  boxsize(-1.0f)
{

}

Particles::~Particles()
{
  delete(particle); particle=0;
}

void Particles::allocate(const index_t n_alloc)
{
  free(particle); np_allocated= 0;
  np_local= np_with_buffers= 0;  id_begin= np_total= 0;
  particle= (Particle*) malloc(sizeof(Particle)*n_alloc);
  assert(particle);
  np_allocated= n_alloc;
}

//
// Halos
//
/*
Halos::Halos() : 
  halos(0), np_allocated(0)
{

}

Halos::~Halos()
{
  free(halos); halos= 0;
}

void Halos::allocate(const int n_alloc)
{
  free(halos);
  halos= (Halos*) malloc(sizeof(Halo)*n_alloc);
  assert(halos);
  np_allocated= n_alloc;
}
*/
