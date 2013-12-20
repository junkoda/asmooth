#include "exchange_buffer.h"
#include <mpi.h>

using namespace std;

inline int get_buffer_index(const Particle p, const float left_threshold, const float right_threshold)
{
  return  (p.x[0] >= left_threshold) + (p.x[0] > right_threshold)
     + 3*((p.x[1] >= left_threshold) + (p.x[1] > right_threshold))
     + 9*((p.x[2] >= left_threshold) + (p.x[2] > right_threshold));
}

void exchange_buffer(MpiInterface const * const mpi, Particles* const particles,
		     const float buffer_factor)
{
  //cout << "np " << particles->np_local << endl; //debug
  //cout << "np_alloc " << particles->np_allocated << endl;

  const float boxsize= particles->boxsize;
  const index_t np= particles->np_local;
  const index_t np_alloc= particles->np_allocated;

  const float left_threshold= buffer_factor*boxsize;
  const float right_threshold= (1.0f-buffer_factor)*boxsize;
  Particle* const p= particles->particle;

  int ibuffer= np_alloc-1;
  //cout << "ibuffer " << ibuffer << endl;
  int nbuf[27];
  for(int i=0; i<27; ++i)
    nbuf[i]= 0;

  //
  // Copy and count buffer particles
  //
  for(int i=0; i<np; ++i) {
    int index= get_buffer_index(p[i], left_threshold, right_threshold);
    if(index != 13) {
      if(ibuffer < np) {
	cerr << "Error: Not enough memory for buffer particles: " 
	     << ibuffer << " " << np << endl;
	abort();
      }
      p[ibuffer--]= p[i];
      nbuf[index]++;
    }
  }

  //debug
  /*
  if(mpi->index() == 0) {
    for(int i=0; i<27; ++i)
      cerr << "dbg: " << nbuf[i] << endl;
  }
  */
 
  //
  // sort buffer particles
  //
  Particle * const q= p + (ibuffer+1); // buffer particles
  int ibuf[27], ibuf_begin[27], ibuf_end[27];
  ibuf[0]= ibuf_begin[0]= 0; ibuf_end[0]= nbuf[0];
  for(int i=1; i<27; ++i) {
    ibuf[i] = ibuf[i-1] + nbuf[i-1];
    ibuf_begin[i]= ibuf[i];
    ibuf_end[i]= ibuf[i] + nbuf[i];
  }  

  const int np_buf= np_alloc - 1 - ibuffer;
  
  // for(int i=0; i<27; ++i)
  //   printf("debug %d %d %d\n", nbuf[i], ibuf[i], ibuf_end[i]);
  // cout << np_buf << " " << ibuf_end[26] << endl;
  assert(np_buf == ibuf_end[26]);
    
  Particle tmp;
  for(int i=0; i<27; ++i) {
    for(int j=ibuf[i]; j<ibuf_end[i]; ++j) {
      int index= get_buffer_index(q[j], left_threshold, right_threshold);
      while(index != i) {
	assert(index > i);
	assert(ibuf[index] < ibuf_end[index]);
	int k= ibuf[index]++;

	assert(0 <= k && k < np_buf);
	tmp= q[k];
	q[k]= q[j];
	q[j]= tmp;

	index= get_buffer_index(tmp, left_threshold, right_threshold);
      }
    }
  }

#ifdef DEBUG
  // debug OK
  mpi->print("debug send buffer\n");
  for(int i=0; i<np_buf; ++i) {
    int index= get_buffer_index(q[i], left_threshold, right_threshold);
    if(!(ibuf_end[index]-nbuf[index] <= i && i < ibuf_end[index]))
      printf("%d %d %d %d\n", i, ibuf_end[index]-nbuf[index],
	     ibuf_end[index], index);
    assert(ibuf_end[index]-nbuf[index] <= i && i < ibuf_end[index]);
  } 
#endif
  
  assert(sizeof(Particle) % sizeof(float) == 0);
  const int particle_len= sizeof(Particle)/sizeof(float);

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

	//printf("%d %d %d\n" ,nbuf[i], (irecv_end-irecv), j);
	//printf("%d %d %d %d %d\n" ,nbuf[i], irecv, irecv_end, (irecv_end-irecv), j);
	//	if(mpi->rank() == 0)
	// cout << "Sendrecv phase " << i << " -> " << j << endl;
	int nrecv= mpi->sendrecv_to(q+ibuf_begin[i], nbuf[i]*particle_len,
				  p+irecv, (irecv_end-irecv)*particle_len, j);
	assert(nrecv % particle_len == 0);
	nrecv /= particle_len;
	assert(irecv + nrecv < irecv_end);

	Particle *b= p+irecv;
	for(int k=0; k<nrecv; ++k) {
	  b->x[0] -= (dn[0]-1)*boxsize;
	  b->x[1] -= (dn[1]-1)*boxsize;
	  b->x[2] -= (dn[2]-1)*boxsize;
	  ++b;
	} // convert coordinate

	irecv += nrecv;
      }
    }
    irecv_end += nbuf[i];
    assert(irecv_end <= np_alloc);
  }

  particles->np_with_buffers= irecv;
  assert(particles->np_with_buffers < particles->np_allocated);

  long np_total= 0;
  long np_long= particles->np_local;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(&np_long, &np_total, 1, MPI_LONG, 
	     MPI_SUM, 0, MPI_COMM_WORLD);
  if(mpi->rank() == 0)
    particles->np_total= np_total;

  //cout << "np_with_buffers " << irecv << endl;
  // np_with_buffers for uniform case is 1.331 np for buffer_factor=0.05
  // buffer_ratio = 6*x + 12*x^2 + 8*x^3
}
