#ifndef EXCHANGE_BUFFER_TEMPLATE_H
#define EXCHANGE_BUFFER_TEMPLATE_H 1

#include <mpi.h>
#include "mpi_interface.h"

using namespace std;

template<class T>
inline int get_buffer_index(const T p, const float left_threshold, const float right_threshold)
{
  return  (p.x[0] >= left_threshold) + (p.x[0] > right_threshold)
     + 3*((p.x[1] >= left_threshold) + (p.x[1] > right_threshold))
     + 9*((p.x[2] >= left_threshold) + (p.x[2] > right_threshold));
}


template<class T>
void exchange_buffer_template(MpiInterface const * const mpi, 
		     ParticleSet<T>* const particles,
		     const float buffer_factor)
{
  const float boxsize= particles->boxsize;
  const float left_threshold= buffer_factor*boxsize;
  const float right_threshold= (1.0f-buffer_factor)*boxsize;
  const index_t np= particles->np_local;
  const index_t np_alloc= particles->np_allocated;
  T* const p= particles->particle;

  index_t ibuffer= np_alloc-1;
  int nbuf[27];
  for(int i=0; i<27; ++i)
    nbuf[i]= 0;

  //
  // Copy and count buffer particles
  //
  for(int i=0; i<np; ++i) {
    int index= get_buffer_index(p[i], left_threshold, right_threshold);
    if(index != 13) {
      assert(ibuffer >= np);
      p[ibuffer--]= p[i];
      nbuf[index]++;
    }
  }

  //
  // sort buffer particles
  //
  T* const q= p + (ibuffer+1); // buffer particles
  int ibuf[27], ibuf_begin[27], ibuf_end[27];
  ibuf[0]= ibuf_begin[0]= 0; ibuf_end[0]= nbuf[0];
  for(int i=1; i<27; ++i) {
    ibuf[i] = ibuf[i-1] + nbuf[i-1];
    ibuf_begin[i]= ibuf[i];
    ibuf_end[i]= ibuf[i] + nbuf[i];
  }  

  const int np_buf= np_alloc - 1 - ibuffer; 
  
  assert(np_buf == ibuf_end[26]);
    
  T tmp;
  for(int i=0; i<27; ++i) {
    for(int j=ibuf[i]; j<ibuf_end[i]; ++j) {
      int index= get_buffer_index(q[j], left_threshold, right_threshold);
      while(index != i) {
	assert(index > i);
	assert(ibuf[index] < ibuf_end[index]);
	int k= ibuf[index]++;

	tmp= q[k];
	q[k]= q[j];
	q[j]= tmp;

	index= get_buffer_index(tmp, left_threshold, right_threshold);
      }
    }
  }

  assert(sizeof(T) % sizeof(float) == 0);
  const int particle_len= sizeof(T)/sizeof(float);

  int irecv= np;
  int irecv_end= ibuffer+1;
  for(int i=0; i<27; ++i) { // buffer source index
    if(i == 13) continue;
    int bn[3];
    bn[2]= i / 9; bn[1]= (i-9*bn[2])/3; bn[0]= i % 3;
    assert(9*bn[2] + 3*bn[1] + bn[0] == i);
    int dn[3];
    for(int j=0; j<27; ++j) { // transfer direction index
      if(j == 13) continue;
      dn[2]= j / 9; dn[1]= (j-9*dn[2])/3; dn[0]= j % 3;
      assert(9*dn[2] + 3*dn[1] + dn[0] == j);

      if(i != 13 && j != 13 &&
	 (dn[0] == 1 || dn[0] == bn[0]) &&
	 (dn[1] == 1 || dn[1] == bn[1]) &&
	 (dn[2] == 1 || dn[2] == bn[2])) {

	int nrecv= mpi->sendrecv_to(q+ibuf_begin[i], nbuf[i]*particle_len,
				  p+irecv, (irecv_end-irecv)*particle_len, j);
	assert(nrecv % particle_len == 0);
	nrecv /= particle_len;

	T *b= p+irecv;
	for(int k=0; k<nrecv; ++k) {
	  b->x[0] -= (dn[0]-1)*boxsize;
	  b->x[1] -= (dn[1]-1)*boxsize;
	  b->x[2] -= (dn[2]-1)*boxsize;
	  ++b;
	}

	irecv += nrecv;
      }
    }
    irecv_end += nbuf[i];
  }

  particles->np_with_buffers= irecv;
}

#endif
