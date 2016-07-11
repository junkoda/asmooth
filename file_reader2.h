#ifndef FILE_READER_H2
#define FILE_READER_H2 1

#include "basic_types.h"
#include "mpi_interface.h"

struct ErrorParticleReader{};

void read_pm_file_xv(const char filebase[], const char redshift[], const int inode, const float buffer_factor, const float shift[], Particles* const particles);

void read_pm_file_zip(const char filebase[], const char redshift[], const int inode, const float buffer_factor, const int nc_node_dim, const int mesh_scale, const float shift[], Particles* const particles);

void write_ascii(const char redshift[], const int inode,  Particles* const particles);

#endif

