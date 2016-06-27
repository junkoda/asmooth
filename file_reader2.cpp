#include <iostream>
#include <cstdio>
#include "basic_types.h"
#include "file_reader.h"

using namespace std;

struct P3m_header
{
  int np_local;
  float a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, 
    cur_projection,cur_halofind,mass_p;
};

struct P3m_zip_header0
{
  int np_local;
  float a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,cur_projection,cur_halofind,mass_p,v_r2i,shake_offset, dummy1 , dummy2;
};


void read_pm_file3(const char filebase[], const char redshift[], const int inode, const float buffer_factor, const int nc_node_dim, const int mesh_scale, const float shift[], Particles* const particles)
{
  //
  // Read PM file and add particles
  //
  // filebase: node or dir/node
  char filename[256];
  Particle* p= particles->particle;
  
  /*
  sprintf(filename, "%s%d/%sxv%d.dat", filebase, inode, redshift, inode);
    
  FILE* const fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to read file: " << filename << endl;
    throw ErrorParticleReader();
  }
  
  //
  // Read header
  //
  P3m_header header;
  const int ret_read_header= fread(&header, sizeof(P3m_header), 1, fp);
  assert(ret_read_header == 1);

  if(particles->np_allocated < header.np_local) {
    cerr << "Not enough space for particles: "
         << particles->np_allocated << " < "
         << header.np_local << endl;
    throw ErrorParticleReader();
  }

  const index_t np_local= header.np_local;
  particles->np_local= np_local;
  particles->np_with_buffers= np_local;
  

  //
  // Read xvfile
  //
  Particle* p= particles->particle;
  const int blocksize=(32*1024*1024)/24;
  const int num_writes=np_local/blocksize+1;

  for(int i=1; i<=num_writes; ++i) {
    const int nplow= (i-1)*blocksize + 1;
    const int nphigh= min(i*blocksize, np_local);
    for(int j=nplow; j<=nphigh; ++j) {
      fread(p[j-1].x, sizeof(float), 3, fp);
      fread(p[j-1].v, sizeof(float), 3, fp);
      p[j-1].x[0] -= shift[0];
      p[j-1].x[1] -= shift[1];
      p[j-1].x[2] -= shift[1];
    }
  }

  const int ret_fclose= fclose(fp);
  assert(ret_fclose == 0);
  */
  
  //
  // Read zip{0123}_<i>.dat
  //
  assert(sizeof(P3m_zip_header0) == 16*sizeof(float));
  
  sprintf(filename, "%s%d/%szip0_%d.dat", filebase, inode, redshift, inode);
  FILE* const fp0= fopen(filename, "r"); assert(fp0);

  sprintf(filename, "%s%d/%szip1_%d.dat", filebase, inode, redshift, inode);
  FILE* const fp1= fopen(filename, "r"); assert(fp1);

  sprintf(filename, "%s%d/%szip2_%d.dat", filebase, inode, redshift, inode);
  FILE* const fp2= fopen(filename, "r"); assert(fp2);

  sprintf(filename, "%s%d/%szip3_%d.dat", filebase, inode, redshift, inode);
  FILE* const fp3= fopen(filename, "r"); assert(fp3);

  P3m_zip_header0 zip_header0;
  float zip_header1[16];
  fread(&zip_header0, sizeof(P3m_zip_header0), 1, fp0);
  fread(zip_header1, sizeof(float), 16, fp1);

  const index_t np_local= zip_header0.np_local;
  particles->np_local= np_local;
  particles->np_with_buffers= np_local;

  
  if(particles->np_allocated < zip_header0.np_local) {
    cerr << "Too many particles to store\n";
    cerr << "np_local= " <<  zip_header0.np_local << endl;
    cerr << "particles->np_allocated= " << particles->np_allocated << endl;
    throw ErrorParticleReader();
  }

  const float v_r2i= zip_header0.v_r2i;
  
  int np_uzip= 0;
  int rr_i4;
  char* rhoc_i1= (char*) &rr_i4; // equivalence(rr_i4, rhoc_i1)
  int xi4[3];
  char * xi1= (char*) xi4;       // equivalence(xi1, xi4);
  short vi2[3];                  // 3 vectors of size 3;
  assert(sizeof(short) == 2);
  assert(sizeof(int) == 4);
  
  for(int k=1; k<=nc_node_dim; ++k) {
    for(int j=1; j<=nc_node_dim; ++j) {
      for(int i=1; i<=nc_node_dim; ++i) {
	rr_i4= 0;

#ifndef BGQ
	// little endian
	fread(rhoc_i1, 1, 1, fp2);     // little endian
#else
	fread(rhoc_i1 + 3, 1, 1, fp2); // big endian
#endif
	if(rr_i4 == 255) { // rr_i4 read throgh rhoc_i1
	  fread(&rr_i4, sizeof(int), 1, fp3);
	}
	for(int l=1; l<=rr_i4; ++l) {
	  xi4[0]= xi4[1]= xi4[2]= 0;
	  np_uzip++;

#ifndef BGQ
	  fread(xi1     , 1, 1, fp0);
	  fread(xi1 +  4, 1, 1, fp0);
	  fread(xi1 +  8, 1, 1, fp0);
#else
	  fread(xi1 +  3, 1, 1, fp0);
	  fread(xi1 +  7, 1, 1, fp0);
	  fread(xi1 + 11, 1, 1, fp0);
#endif
	  fread(vi2, sizeof(short), 3, fp1);
	  int index= np_uzip - 1;
	  assert(0 <= index && index < np_local);
	  p[index].x[0]=
	    mesh_scale*((xi4[0] + 0.5f)/256.0f + i - 1.0f) - shift[0];
	  p[index].x[1]=
	    mesh_scale*((xi4[1] + 0.5f)/256.0f + j - 1.0f) - shift[1];
	  p[index].x[2]=
	    mesh_scale*((xi4[2] + 0.5f)/256.0f + k - 1.0f) - shift[2];
	  
	  p[index].v[0] = vi2[0] / v_r2i;
	  p[index].v[1] = vi2[1] / v_r2i;
	  p[index].v[2] = vi2[2] / v_r2i;
	}
      }
    }
  }

  // Check end of file)
  char test_i1;
  assert(feof(fp0) == 0);     // should not be EOF yet
  fread(&test_i1, 1, 1, fp0);
  assert(feof(fp0) != 0);     // should be EOF

  assert(feof(fp1) == 0);
  fread(&test_i1, 1, 1, fp1);
  assert(feof(fp1) != 0);

  printf("np_uzip = %d\n", np_uzip);
  printf("nc_node_dim = %d\n", nc_node_dim);
  
  if(np_uzip != np_local) {
    cerr << "Something wrong with reading dark matter zipped files\n";
    throw ErrorParticleReader();
  }

  particles->boxsize= mesh_scale*nc_node_dim;
  
  fclose(fp0);
  fclose(fp1);
  fclose(fp2);
}

