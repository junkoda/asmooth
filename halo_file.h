#ifndef HALO_FILE_H
#define HALO_FILE_H 1

#include "basic_types.h"
#include "mpi_interface.h"

/*
int read_halo_file(const char filename[], ParticleSet<Halo>* const halos,
		   const float shift[]);
*/

void read_and_exchage_halo(const char filename[],
			   const float z, const float omega_m,
			   MpiInterface const * const mpi,
			   ParticleSet<Halo>* const halos,
			   const float buffer_factor,
			   const float shift[]);
/*
int read_halo_file_all(const char redshift[], 
		       const int file_index, const int file_nc,
		       ParticleSet<Halo>* const halos,
		       const float buffer_factor,
		       const float shift[]);
*/
#endif
