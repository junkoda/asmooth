#ifndef COARSE_MESH_H
#define COARSE_MESH_H 1

void clear_mesh(float mesh[], const int nc);
double mesh_total(const float mesh[], const int nc);
void assign_on_mesh(Particles const * const particles,
		    float mesh[], const int nc,
		    const int sign= 1);
void write_mesh(char filename[], float* mesh, const int nc, const float boxsize);

void write_mesh_separate(const char base_filename[], const int inode, 
			 float *const mesh, const int nc, const float boxsize);

void assign_on_density_mesh(Particles const * const particles,
		    float mesh[], const int nc,
		    const int sign= 1);
#endif
