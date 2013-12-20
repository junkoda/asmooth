#ifndef NBR_SEARCH_H
#define NBR_SEARCH_H 1

#include "basic_types.h"
#include "kdtree_balanced.h"
#include "k_neighbors.h"

void set_nbr_search(Node* const root_, Particle* const p0);
void neighbor_search(Particle* const p, KNeighbors* const knbrs);

int mark_halo_particles(Node const * const n, const float x[], const float r);
#endif
