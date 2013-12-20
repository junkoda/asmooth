//
// Neighbor Search Rountines: find fixed number of nbrs and calc density
//

#include <iostream>
#include <cmath>
#include <assert.h>
#include "nbr_search.h"
#include "k_neighbors.h"

using namespace std;

static Node * root;
static Particle * particle;

void find_nbr_recursive(Node const * const,
		       const float x[], KNeighbors* const knbrs);

inline float norm2(const float x[], const float y[]);
float norm2(const float x[], const float y[]) {
  const float dx= x[0]-y[0];
  const float dy= x[1]-y[1];
  const float dz= x[2]-y[2];
  return dx*dx + dy*dy + dz*dz;
}

inline float ker(const float x[], const float y[], const float hinv)
{
  //debug!
  assert(abs(x[0]-y[0])*hinv <= 1.0f);
  assert(abs(x[1]-y[1])*hinv <= 1.0f);
  assert(abs(x[2]-y[2])*hinv <= 1.0f);

  return hinv*hinv*hinv*(1.0f-abs(x[0]-y[0])*hinv)
                       *(1.0f-abs(x[1]-y[1])*hinv)
                       *(1.0f-abs(x[2]-y[2])*hinv);
}

inline float spline_ker(const float x[], const float y[], const float hinv)
{
  const float fac= 8.0f/M_PI*hinv*hinv*hinv;
  const float r= sqrt(norm2(x,y))*hinv;
  //if(r < 0.0f || r > 1.0f)
  //  cerr << "r= " << r << endl;
  assert(r >= 0.0f && r <= 1.00001f);
  if(r < 0.5f)
    return fac*(1 - 6*r*r + 6*r*r*r);
  //else if(r < 1.0f) {
  float r1= 1.0f-r;
  return fac*(2*r1*r1*r1);
    //}
    //return 0.0f;
}

void set_nbr_search(Node* const root_, Particle* const p0)
{
  root= root_;  // root node
  particle= p0; // 0th particle
}

void neighbor_search(Particle* const p, KNeighbors* const knbrs)
{
  assert(root); assert(particle);
  knbrs->clear();
  //float const * const x= p->x;
  find_nbr_recursive(root, p->x, knbrs);

  // Assign density at particle positions
  if(knbrs->get_kth_value() <= 0.0f) {
    cerr << knbrs->get_kth_value() << endl;
    assert(false);
  }
  p->rk= sqrt(knbrs->get_kth_value());
  assert(p->rk > 0.0f);
  const float hinv= 1.0f/p->rk;
  const int k= knbrs->k;

  for(int i=1; i<k; ++i) {
    Particle* const q= knbrs->nbr(i);
    //if(q->id != p->id) {
    q->dens += ker(p->x, q->x, hinv);      // tsc
      //q->dens += spline_ker(p->x, q->x, hinv); // spline kernel
      //}
  }


  //cout << knbrs->get_kth_value() << endl;
  //return sqrt(knbrs->get_kth_value());
}

void find_nbr_recursive(Node const * const n, 
			const float x[], KNeighbors* const knbrs)
{
  const float r= sqrt(knbrs->get_kth_value());
  const float y= x[n->direction];

  if(y+r < n->left || y-r > n->right)
    return;

  if(n->content[0] < 0) {
    for(index_t i= -n->content[0]-1; i<-n->content[1]; ++i) {
      float d2= norm2(x, particle[i].x);
      if(d2 < r*r)
	knbrs->push(d2, particle+i);
    }
    return;
  }
  
  const int subtree_cases= (n->content[0] > 0)
                           + 2*(n->content[1] > 0);
  switch(subtree_cases) {
  case 0:
    assert(false); break;
  case 1:
    find_nbr_recursive(root + n->content[0], x, knbrs); break;
  case 2:
    find_nbr_recursive(root + n->content[1], x, knbrs); break;
  case 3:
    Node const * const next_node= root + n->content[0];
    const int near_subtree= x[next_node->direction] >= next_node->right;
    find_nbr_recursive(root + n->content[near_subtree],
		       x, knbrs);
    find_nbr_recursive(root + n->content[!near_subtree],
		       x, knbrs);
    break;
  }
}
  
//
//
//
/*
void mark_halo_particles(const float x[], const float r)
{
  assert(root); assert(particle);
  mark_halo_recursive(root, x, r);
}
*/

int mark_halo_particles(Node const * const n, const float x[], const float r)
{
  int np_halo= 0;
  const float y= x[n->direction];

  if(y+r < n->left || y-r > n->right)
    return np_halo;

  if(n->content[0] < 0) {
    for(index_t i= -n->content[0]-1; i<-n->content[1]; ++i) {
      float d2= norm2(x, particle[i].x);
      if(d2 < r*r) {
	particle[i].igrp= -1;
	np_halo++;
      }
    }
    return np_halo;
  }
  
  if(n->content[0])
    np_halo += mark_halo_particles(root + n->content[0], x, r);
  if(n->content[1])
    np_halo += mark_halo_particles(root + n->content[1], x, r);

  return np_halo;
}
