#ifndef FILE_READER_H2
#define FILE_READER_H2 1

#include "basic_types.h"
#include "mpi_interface.h"

struct ErrorParticleReader{};

void read_pm_file3(const char filebase[], const char redshift[], const int inode, const float buffer_factor, const int nc_node_dim, const int mesh_scale, const float shift[], Particles* const particles);

#endif

