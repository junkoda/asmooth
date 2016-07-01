#ifndef KDTREE_BALANCED_H
#define KDTREE_BALANCED_H 1

#include "basic_types.h"

struct Node {
  int direction;
  float left, right;
  index_t content[2];
};

class KDTree {
 public:
  KDTree();
  virtual ~KDTree();
  virtual void allocate(const int np_max, const int quota= 16)=0;
  virtual index_t construct(Particles* const, const float, const float)=0;
  size_t estimate_mem(const int np_max, const int quota_);
  Node* root() const { return my_nodes; }
  Particles* particles(){ return my_particles; }
 protected:
  Node* my_nodes;
  Particles* my_particles;
  int nn_alloc; // number of nodes allocated
  int quota;
};

class KDTreeBalanced : public KDTree {
 public:
  virtual void allocate(const int np_max, const int quota= 16);
  virtual index_t construct(Particles* const, const float, const float);  
 private:
  index_t split_particles(const int inode, 
			  const float left[], const float right[],
			  const index_t ibegin, const index_t iend,
			  const int quota, const int d);
};

class KDTreeSimple : public KDTree {
 public:
  virtual void allocate(const int np_max, const int quota= 16);
  virtual index_t construct(Particles* const, const float, const float);  
 private:
  index_t construct_tree_recursive(index_t inode,
				   const index_t particle_begin,
				   const index_t particle_end,
				   const int this_direction,
				   const float left3[],
				   const float right3[]);
};



#endif
