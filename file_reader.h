#ifndef FILE_READER_H
#define FILE_READER_H 1

#include "basic_types.h"
#include "mpi_interface.h"

struct ErrorParticleReader{};

//void read_pm_file(MpiInterface const * const mpi, const char filename[],
//		  Particles* const);
void read_pm_file(const char filename[], Particles* const);
void read_pm_file2(const char filename[], Particles* const, 
		   const float shift[]);

void read_pm_file_all(const char redshift[], 
		      const int file_index, const int file_nc,
		      Particles* const particles,
		      const float buffer_factor,
		      const float shift[]);

//void read_shift(const char redsfhit[], float shift[]);

#endif

