#ifndef FOF_H
#define FOF_H 1

#include "basic_types.h"
#include "kdtree_balanced.h"

index_t group_fof(Particles* const particles, 
	       KDTree const * const tree,
	       const float linking_length /* linking length */);

index_t remove_small_groups(Particles* const particles,
			    int* const nfriends, const int ngrp);

void write_fof_halos(const char filename[],
		     Particles* const particles,
		     void* const buf, const int nhalo);
#endif
