#ifndef EXCHANGE_BUFFER_H
#define EXCHANGE_BUFFER_H 1

#include "basic_types.h"
#include "mpi_interface.h"

void exchange_buffer(MpiInterface const * const, Particles * const,
		     const float buffer_factor);


#endif
