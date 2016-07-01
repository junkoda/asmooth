#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>
#include "kdtree_balanced.h"

using namespace std;

/*
struct Box {
  float left[3], right[3];
};
*/

//
// KDTree
//
KDTree::KDTree() :
  my_nodes(0), my_particles(0), nn_alloc(0), quota(0)
{

}

KDTree::~KDTree()
{
  free(my_nodes);
}

size_t KDTree::estimate_mem(const int np_max, const int quota_)
{
  index_t n_leaf= 1;
  
  while(16*n_leaf < np_max)
    n_leaf <<= 1;

  nn_alloc= 2*n_leaf - 1;

  return sizeof(Node)*(size_t)(nn_alloc);
}



//
// KDTreeBalanced
//  
void KDTreeBalanced::allocate(const int np_max, const int quota_)
{
  index_t n_leaf= 1;
  quota= quota_;
  assert(quota>0);
  
  while(16*n_leaf < np_max)
    n_leaf <<= 1;
  //cout << "n_leaf " << np_max << " " <<  n_leaf << endl;

  nn_alloc= 2*n_leaf - 1;

  my_nodes= (Node*) malloc(sizeof(Node)*(size_t)(nn_alloc));
  assert(my_nodes);
}

index_t KDTreeBalanced::construct(Particles* const particles,
			       const float left, const float right)
{
  assert(my_nodes); assert(nn_alloc);
  my_particles= particles;

  //Box box;
  //box.left[0]= box.left[1]= box.left[2]= left;
  //box.right[0]= box.right[1]= box.right[2]= right;

  //int inode= 0;
  int ibegin= 0, iend= particles->np_with_buffers;
  float left3[]= {left, left, left};
  float right3[]= {right, right, right};

  split_particles(0, left3, right3, ibegin, iend, quota, 0);
  return nn_alloc;
}

index_t KDTreeBalanced::split_particles(const int inode, 
				    const float left3[], const float right3[],
				    const index_t ibegin, const index_t iend,
				    const int quota, const int this_direction)
{
  Particle* const p= my_particles->particle;
  Particle tmp;
  const int i= (this_direction+1)%3; // cutting direction

  //float xl= left3[i], xr= right3[i];
  index_t ib= ibegin, ie= iend;

  index_t ic= ibegin+(iend-ibegin)/2;

  Node* const n= my_nodes+inode;
  n->direction= this_direction;
  n->left= left3[this_direction]; n->right= right3[this_direction];

  if(iend - ibegin <= quota) {
    n->content[0]= -ibegin-1;
    n->content[1]= -iend;
    return inode;
  }

  //n->content[0]= inode << 1 + 1;
  //n->content[1]= inode << 2 + 2;

  float xc;
  while(1) {

    //if(++ithloop < 10)
    //  xc= xl + 0.5f*(xr-xl);
    //else
      xc= p[ib].x[i];

    //xc= p[ib].x[i];
    //cout << xc << endl;
    assert(ib < ie);

    int ileft=ib, iright=ie-1;
    while(1) {
      while(ileft<ie && p[ileft].x[i] < xc) ++ileft;
      while(iright>=ib && p[iright].x[i] >= xc) --iright;
      if(ileft < iright) {
	tmp= p[ileft]; p[ileft]= p[iright]; p[iright]= tmp; // swap
	ileft++; iright--;
      }
      else
	break;
    }
    //printf("rl %d %d %d %d %d %f\n", ileft, iright, ic, ib, ie, xc);

    if(iright == ib-1) {  //&& ithloop >= 10) {
      if(ib == ic)
	break;
      ib++;
      continue;
    }

    else if(ileft-iright != 1) {
      for(int j=iright; j<=ileft; ++j)
	printf("%8d %13.8e\n", j, p[j].x[i]);
    }
    assert(ileft-iright == 1); // true only if all particles have different x

    const int ic_tmp= ileft;
    if(ic_tmp < ic) {
      ib= ic_tmp; //xl= xc;
    }
    else if(ic_tmp > ic) {
      ie= ic_tmp; //xr= xc;
    }
    else
      break;
    
  }

  //debug: expensive test
  //assert(xc == p[ic].x[i]);

  for(int j= ibegin; j<ic; ++j) {
    if(!(p[j].x[i] <= xc))
      printf("#1 %d %e %e\n", j, p[j].x[i], xc);
    assert(p[j].x[i] <= xc);
  }
  for(int j=ic; j<iend; ++j) {
    if(!(p[j].x[i] >= xc))
      printf("#2 %d %13.8e %13.8e %13.8e\n", j, p[j].x[i], p[ic].x[i], xc);
    assert(p[j].x[i] >= xc);
  }

  //cout << "ic = " << ic << endl;
  float left3_next[]= {left3[0], left3[1], left3[2]};
  float right3_next[]= {right3[0], right3[1], right3[2]};
  left3_next[i]= xc;
  right3_next[i]= xc;


  n->content[0]= inode+1;
  index_t node_l_end= 
    split_particles(n->content[0], left3, right3_next, ibegin, ic, quota, i);
  n->content[1]= node_l_end+1;
  index_t node_r_end=
    split_particles(n->content[1], left3_next, right3, ic, iend, quota, i);


    
  //else
  //  break;
  return node_r_end;
}

//
// KDTreeSimple (kdtree bisecting space in half
//  

void KDTreeSimple::allocate(const int np_max, const int quota_)
{
  index_t n_leaf= 1;
  quota= quota_;
  assert(quota>0);
  
  while(16*n_leaf < np_max)
    n_leaf <<= 1;
  //cout << "n_leaf " << np_max << " " <<  n_leaf << endl;

  nn_alloc= 2*n_leaf - 1;
  my_nodes= (Node*) malloc(sizeof(Node)*(size_t)(nn_alloc)*2);
  assert(my_nodes);
}

index_t KDTreeSimple::construct(Particles* const particles,
			     const float left, const float right)
{
  assert(my_nodes); assert(nn_alloc);
  my_particles= particles;

  //Box box;
  //box.left[0]= box.left[1]= box.left[2]= left;
  //box.right[0]= box.right[1]= box.right[2]= right;
  float left3[]= {left, left, left};
  float right3[]= {right, right, right};

  //int inode= 0;
  int ibegin= 0, iend= particles->np_with_buffers;

  return construct_tree_recursive(0, ibegin, iend, 0, left3, right3);
}

index_t KDTreeSimple::construct_tree_recursive(index_t inode,
					       const index_t particle_begin,
					       const index_t particle_end,
					       const int this_direction,
					       const float left3[],
					       const float right3[]) 
{
  assert(particle_begin < particle_end);

  Particle* const p= my_particles->particle;
  Node* n= my_nodes + inode;
  n->direction= this_direction;
  n->left= left3[this_direction];
  n->right= right3[this_direction];

  if(particle_end - particle_begin < quota) {
    n->content[0]= -particle_begin-1;
    n->content[1]= -particle_end;
    return inode;
  }

  n->content[0]= n->content[1]= 0;
  const int i= (this_direction + 1) % 3; // cutting direction
  const float mid= left3[i]+0.5f*(right3[i]-left3[i]); 

  // split_particles
  int j1= particle_begin;
  int j2= particle_end-1;
  
  while(1) {
    while(p[j1].x[i] < mid && j1 < particle_end)
      ++j1;
    while(p[j2].x[i] >= mid && j2 >= particle_begin)
      --j2;
    if(j1 < j2) {
      const Particle temp= p[j1];
      p[j1]= p[j2];
      p[j2]= temp;
      ++j1; --j2;
    }
    else
      break;
  }
  assert(j1 >= particle_begin && j1 <= particle_end && (j1 - j2 == 1));

  const index_t particle_mid= j1;
  
  if(particle_mid != particle_begin) {
    float right_next[]= {right3[0], right3[1], right3[2]};
    right_next[i]= mid;
    n->content[0]= inode+1;
    inode= construct_tree_recursive(n->content[0],
				    particle_begin, particle_mid,
				    i, left3, right_next);
  }
  if(particle_mid < particle_end) {
    float left_next[]= {left3[0], left3[1], left3[2]};
    left_next[i]= mid;
    n->content[1]= inode+1;
    inode= construct_tree_recursive(n->content[1],
				    particle_mid, particle_end,
				    i%3, left_next, right3);
  }
  
  assert(n->content[0] > 0 || n->content[1] > 0);
  return inode;
}

