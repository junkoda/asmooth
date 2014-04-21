#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>

using namespace std;


int main(int argc, char* argv[])
{
  assert(sizeof(int) == 4);

  if(argc != 5) {
    cout << "reformat <nc_cpu> <out_dir> <redshift> <in_dir>\n"
         << "\treformat 8 7.221\n";
    return 0;
  }

  //const int nc= atoi(argv[1]);
  const int nc_cpu= atoi(argv[1]);

  //float* buf= 0;//= (float*) malloc(sizeof(float)*nc);
  float* buf_n= 0;
  float* buf_h= 0;
  float* buf_v= 0;

  //cout << argv[argc-1] << " " << nc_cpu << endl;

  //
  // Open large file
  //
  char filename_n[128];
  char filename_h[128];
  char filename_v[128];
  //char const * const file_pattern= argv[2];
  //if(argc == 3)
  //  sprintf(filename, "%s_all.dat", argv[argc-1]);
  //else
  sprintf(filename_n, "%s/%sn_all.dat", argv[2], argv[3]);
  sprintf(filename_h, "%s/%sntot_all.dat", argv[2], argv[3]);
  sprintf(filename_v, "%s/%sv_all.dat", argv[2], argv[3]);

  FILE* fp_n_out= fopen(filename_n, "w"); assert(fp_n_out);
  FILE* fp_h_out= fopen(filename_h, "w"); assert(fp_h_out);
  FILE* fp_v_out= fopen(filename_v, "w"); assert(fp_v_out);

  cout << "reformat2 " << filename_n << endl;

  int count[nc_cpu*nc_cpu];
  long double total1= 0.0, total2= 0.0;
  
  for(int i=0; i<nc_cpu*nc_cpu; ++i)
    count[i]=0;

  int ret;

  for(int iz_cpu=0; iz_cpu<nc_cpu; ++iz_cpu) {
    //
    // For each xy slice of cpu plane
    //
    int index_cpu_begin= nc_cpu*nc_cpu*iz_cpu;
    FILE *fp_n[nc_cpu*nc_cpu], *fp_h[nc_cpu*nc_cpu], *fp_v[nc_cpu*nc_cpu];;

    //
    // Open small files
    //
    int nc=0, ncd;
    for(int index_cpu= 0; index_cpu<nc_cpu*nc_cpu; ++index_cpu) {
      sprintf(filename_n, "%s/%sn%d.dat", argv[4], argv[3], index_cpu_begin+index_cpu);
      fp_n[index_cpu] = fopen(filename_n, "r"); // assert(fp_n[index_cpu]);
      if(fp_n[index_cpu] == 0) {
	cerr << "Unable to open: " << filename_n << endl; abort();
      }

      sprintf(filename_h, "%s/%sh%d.dat", argv[4], argv[3], index_cpu_begin+index_cpu);
      fp_h[index_cpu] = fopen(filename_h, "r"); //assert(fp_h[index_cpu]);
      if(fp_h[index_cpu] == 0) {
	cerr << "Unable to open: " << filename_h << endl; abort();
      }

      sprintf(filename_v, "%s/%svel%d.dat", argv[4], argv[3], index_cpu_begin+index_cpu);
      fp_v[index_cpu] = fopen(filename_v, "r"); //assert(fp_v[index_cpu]);
      if(fp_v[index_cpu] == 0) {
	cerr << "Unable to open: " << filename_v << endl; abort();
      }
      
      fread(&nc, sizeof(int), 1, fp_n[index_cpu]);
      assert(nc > 0);
      //cerr << "nc " << nc << endl;

      fread(&ncd, sizeof(int), 1, fp_h[index_cpu]);
      assert(ncd == nc);

      fread(&ncd, sizeof(int), 1, fp_v[index_cpu]);
      assert(ncd == nc);
    }

    if(buf_n == 0) {
      buf_n= (float*) malloc(sizeof(float)*nc); assert(buf_n);
      buf_h= (float*) malloc(sizeof(float)*nc); assert(buf_h);
      buf_v= (float*) malloc(sizeof(float)*nc*3); assert(buf_v);
    }

    if(iz_cpu == 0) {
      int nc_global= nc*nc_cpu;
      int nc_header[]= {nc_global, nc_global, nc_global};
      fwrite(nc_header, sizeof(int), 3, fp_n_out);
      fwrite(nc_header, sizeof(int), 3, fp_h_out);
      fwrite(nc_header, sizeof(int), 3, fp_v_out);
    }

    for(int iz=0; iz<nc; ++iz) {
      //cerr << "iz " << iz << endl;
      for(int iy_cpu=0; iy_cpu<nc_cpu; ++iy_cpu) {
	for(int iy=0; iy<nc; ++iy) {
	  for(int ix_cpu=0; ix_cpu<nc_cpu; ++ix_cpu) {
	    int index_cpu= iy_cpu*nc_cpu + ix_cpu;
	    count[index_cpu]++;
	    ret= fread(buf_n, sizeof(float), nc, fp_n[index_cpu]);
	    assert(ret == nc);
	    ret= fread(buf_h, sizeof(float), nc, fp_h[index_cpu]);
	    assert(ret == nc);
	    ret= fread(buf_v, sizeof(float), nc*3, fp_v[index_cpu]);
	    assert(ret == nc*3);
	    
	    for(int j=0; j<nc; ++j) {
	      float n= buf_n[j]/8.0f;
	      total1 += n;
	      buf_h[j] += buf_n[j];
	      total2 += buf_h[j]/8.0f;
	      buf_v[3*j  ] *= n;
	      buf_v[3*j+1] *= n;
	      buf_v[3*j+2] *= n;
	    }

	    ret= fwrite(buf_n, sizeof(float), nc, fp_n_out);
	    assert(ret == nc);

	    ret= fwrite(buf_h, sizeof(float), nc, fp_h_out);
	    assert(ret == nc);

	    ret= fwrite(buf_v, sizeof(float), nc*3, fp_v_out);
	    assert(ret == nc*3);

	    
	    //for(int i=0; i<nc; ++i)
	    //  total += buf[i];
	  }
	}
      }
    }

    for(int index_cpu= 0; index_cpu<nc_cpu*nc_cpu; ++index_cpu) {
      //assert(feof(fp[index_cpu]));
      assert(count[index_cpu] == nc*nc);
      count[index_cpu]= 0;
      ret= fclose(fp_n[index_cpu]); assert(ret == 0);
      ret= fclose(fp_h[index_cpu]); assert(ret == 0);
      ret= fclose(fp_v[index_cpu]); assert(ret == 0);
    }
  }

  ret= fclose(fp_n_out); assert(ret == 0);
  ret= fclose(fp_h_out); assert(ret == 0);
  ret= fclose(fp_v_out); assert(ret == 0);

  printf("total %s %15.4Lf %15.4Lf\n", argv[3], total1, total2);
  //cout << "total " << total1 << " " << total2 << endl;

  return 0;
}
