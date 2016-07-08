#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>

using namespace std;


int main(int argc, char* argv[])
{
  assert(sizeof(int) == 4);
  if(argc !=5) {
    cout << "reformat <nc_cpu> <out_dir> <file_base> <n_components>\n"
         << "\treformat 8 7.221c 1\n"
         << "\treformat 8 7.221v 3\n";
    return 0;
  }

  //const int nc= atoi(argv[1]);
  const int nc_cpu= atoi(argv[1]);
  const int n_component= atoi(argv[4]);

  float* buf= 0;//= (float*) malloc(sizeof(float)*nc);

  //cout << argv[argc-1] << " " << nc_cpu << endl;

  //
  // Open large file
  //
  char filename[128];
  //char const * const file_pattern= argv[2];
  //if(argc == 3)
  //  sprintf(filename, "%s_all.dat", argv[argc-1]);
  //else
  sprintf(filename, "%s/%s_all.dat", argv[2], argv[3]);

  FILE* fp_out= fopen(filename, "w");
  if(fp_out == 0) {
    cerr << "Unable to write to file: " << filename << endl;
    return 1;
  }
  cout << filename << endl;

  int count[nc_cpu*nc_cpu];
  long double total= 0.0;
  
  for(int i=0; i<nc_cpu*nc_cpu; ++i)
    count[i]=0;
  for(int iz_cpu=0; iz_cpu<nc_cpu; ++iz_cpu) {
    //
    // For each xy slice of cpu plane
    //
    int index_cpu_begin= nc_cpu*nc_cpu*iz_cpu;
    FILE* fp[nc_cpu*nc_cpu];

    //
    // Open small files
    //
    int nc=0;
    for(int index_cpu= 0; index_cpu<nc_cpu*nc_cpu; ++index_cpu) {
      sprintf(filename, "%s%d.dat", argv[3], index_cpu_begin+index_cpu);
      fp[index_cpu] = fopen(filename, "r");
      if(fp[index_cpu] == 0) {
	cerr << "Unable to open: " << filename << endl;
	return 1;
      }
      fread(&nc, sizeof(int), 1, fp[index_cpu]);
      assert(nc > 0);
    }

    const int ndata= nc*n_component;

    if(buf == 0) {
      buf= (float*) malloc(sizeof(float)*ndata);
      if(buf == 0) {
	cerr << "Unable to allocate memory\n";
	return 1;
      }
    }

    if(iz_cpu == 0) {
      int nc_global= nc*nc_cpu;
      int nc_header[]= {nc_global, nc_global, nc_global};
      fwrite(nc_header, sizeof(int), 3, fp_out);
    }

    for(int iz=0; iz<nc; ++iz) {
      for(int iy_cpu=0; iy_cpu<nc_cpu; ++iy_cpu) {
	for(int iy=0; iy<nc; ++iy) {
	  for(int ix_cpu=0; ix_cpu<nc_cpu; ++ix_cpu) {
	    int index_cpu= iy_cpu*nc_cpu + ix_cpu;
	    count[index_cpu]++;
	    if(fread(buf, sizeof(float), ndata, fp[index_cpu]) != ndata) {
	      cerr << "Unable to read data: " << 
		iz_cpu << " " << index_cpu << "\n";
	      return 1;
	    }
	    if(fwrite(buf, sizeof(float), ndata, fp_out) != ndata) {
	      cerr << "Fail to write large data\n";
	      return 1;
	    }
	    for(int i=0; i<nc; ++i)
	      total += buf[i];
	  }
	}
      }
    }

    for(int index_cpu= 0; index_cpu<nc_cpu*nc_cpu; ++index_cpu) {
      //assert(feof(fp[index_cpu]));
      assert(count[index_cpu] == nc*nc);
      count[index_cpu]= 0;
      int ret= fclose(fp[index_cpu]);
      assert(ret == 0);
    }
  }

  fclose(fp_out);
  cout << "Refomat finished successfully. sizeof(int) = " << sizeof(int) << endl;
  cout << "total " << total << endl;

  return 0;
}
