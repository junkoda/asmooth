#include <iostream>
#include <cmath>
#include <queue>
#include <vector>
#include "fof.h"
#include "kdtree_balanced.h"

using namespace std;

static float ll;
static Particle* particle;
static Node* root;

inline float norm2(const float x[], const float y[]) {
  const float dx= x[0]-y[0];
  const float dy= x[1]-y[1];
  const float dz= x[2]-y[2];
  return dx*dx + dy*dy + dz*dz;
}

struct FofHalo {
  float x[3], v[3];
  int nfof;
};

void link_friends_recursive1(Node const * const n,
			     const float x[],
			     queue<index_t>* const friends,
			     const index_t igrp);

index_t check_tree_recursive(Node* const n, index_t np);

index_t group_fof(Particles* const particles, 
	       KDTree const * const tree,
	       const float linking_length /* linking length */)
{
  root= tree->root();
  particle= particles->particle;
  ll= linking_length;

  index_t ngroup= 0;

  queue<index_t> friends;

#ifdef DEBUG
  const index_t np_node= check_tree_recursive(root, 0);
  cout << "debug: check_tree_recursive: " << np_node << "\n";
#endif

  //int i=0; // for all particles
  //if(particle[i].igrp > 0) continue;
  const int np= particles->np_with_buffers;
  for(int i=0; i<np; ++i) {
    if(particle[i].igrp > 0)
      continue;
    ngroup++;
    particle[i].igrp= ngroup;
    friends.push(i);
    while(!friends.empty()) {
      float* const x= particle[friends.front()].x;
      friends.pop();
      link_friends_recursive1(tree->root(),
			      x,
			      &friends,
			      ngroup);
    }
  }
   
  return ngroup;
}

void link_friends_recursive1(Node const * const n,
			    const float x[],
			    queue<index_t>* const friends,
			    const index_t igrp)
{
  const float y= x[n->direction];

  if(y+ll < n->left || y-ll > n->right)
    return;

  if(n->content[0] < 0) {
    for(index_t i= -n->content[0]-1; i<-n->content[1]; ++i) {
      if(particle[i].igrp < igrp) {
	float d2= norm2(x, particle[i].x);
	if(d2 <= ll*ll) {
	  assert(particle[i].igrp == 0);
	  friends->push(i);
	  particle[i].igrp= igrp;
	}
      }
    }
    return;
  }
  
  if(n->content[0])
    link_friends_recursive1(root + n->content[0], x, friends, igrp);
  if(n->content[1])
    link_friends_recursive1(root + n->content[1], x, friends, igrp);

}

index_t remove_small_groups(Particles* const particles,
			    int* const nfriends, const int ngrp)
{
  const int threshold= 20;

  for(int i=0; i<=ngrp; ++i)
    nfriends[i]= 0;

  const int np= particles->np_with_buffers;
  Particle * const p= particles->particle;
  assert(ngrp <= np);

  for(int j=0; j<np; ++j) {
    //assert(p[j].igrp < ngrp); // trivial
    nfriends[p[j].igrp]++;
    //int i= p[j].igrp;
    //nfriends[i]++;
  }

  /*
  for(int j=0; j<np; ++j) {
    if(nfriends[p[j].igrp] < threshold)
      p[j].igrp= 0;
  }
  */

  // remap [1,...,ngrp] -> [0,...,nhalo]
  index_t nhalo= 0;
  int * const halo_index= nfriends;

  for(int j=1; j<=ngrp; ++j) {
    if(nfriends[j] >= threshold)
      halo_index[j]= ++nhalo;
    else
      halo_index[j]= 0;
  }
  assert(nfriends[0] == 0);

  for(int j=0; j<np; ++j) {
    //if(nfriends[p[j].igrp] < threshold)
    p[j].igrp= nfriends[p[j].igrp];
  }

  /*
#ifdef DEBUG
  vector<index_t> v;
  v.reserve(nhalo);
  for(int j=1; j<=ngrp; ++j) {
    if(nfriends[j] >= threshold)
      v.push_back(nfriends[j]);
  }
  sort(v.begin(), v.end());
  FILE* fp= fopen("friends.txt", "w");
  assert(fp);
  for(int j=0; j<v.size(); ++j)
    fprintf(fp, "%d %d\n", j, v[v.size()-j-1]);
  fclose(fp);
  // debug end
#endif
  */

  return nhalo;
}


void write_fof_halos(const char filename[],
			Particles* const particles,
			void* const buf, const int nhalo)
{
  const int np= particles->np_with_buffers;
  assert(sizeof(FofHalo)/sizeof(float)*(nhalo) < np); // sizeof buffer is np

  FofHalo* const h= (FofHalo*) buf;
  for(int i=0; i<=nhalo; ++i) {
    h[i].x[0]= h[i].x[1]= h[i].x[2]= 0.0f;
    h[i].v[0]= h[i].v[1]= h[i].x[2]= 0.0f;
    h[i].nfof= 0;
  }

  Particle * const p= particles->particle;

  for(int j=0; j<np; ++j) {
    int i= p[j].igrp - 1;
    if(i >= 0) {
      h[i].nfof ++;
      for(int k=0; k<3; ++k) {
	h[i].x[k] += p[j].x[k];
	h[i].v[k] += p[j].v[k];
      }
    }
  }

  int nlocal_halo= 0;
  const float boxsize= particles->boxsize;
  for(int j=0; j<nhalo; ++j) {
    for(int k=0; k<3; ++k) {
      h[j].x[k] /= h[j].nfof;
      h[j].v[k] /= h[j].nfof;
    }
    if(0.0f < h[j].x[0] && h[j].x[0] < boxsize &&
       0.0f < h[j].x[1] && h[j].x[1] < boxsize &&
       0.0f < h[j].x[2] && h[j].x[2] < boxsize) {
      h[nlocal_halo++]= h[j];
    }
  }

  FILE* fp= fopen(filename, "w");
  assert(fp);
  fwrite(&nlocal_halo, sizeof(int), 1, fp);
  int size_halo= sizeof(FofHalo);
  fwrite(&size_halo, sizeof(int), 1, fp);
  //fwrite(h, sizeof(FofHalo), nhalo, fp); // wrong line, used for 425Mpc, 1Gpc
  fwrite(h, sizeof(FofHalo), nlocal_halo, fp);
  fwrite(&nlocal_halo, sizeof(int), 1, fp);
  fclose(fp);
}

index_t check_tree_recursive(Node* const n, index_t np)
{
  // debug
  //const float y= x[n->direction];

  if(n->content[0] < 0) {
    for(index_t i= -n->content[0]-1; i<-n->content[1]; ++i) {
      assert(n->left <= particle[i].x[n->direction]);
      assert(n->right >= particle[i].x[n->direction]);
      np++;
    }
    return np;
  }
  
  if(n->content[0])
    np = check_tree_recursive(root + n->content[0], np);
  if(n->content[1])
    np = check_tree_recursive(root + n->content[1], np);

  return np;
}
