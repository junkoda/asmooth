#include "open.h"

using namespace std;

FILE* open(const char filename[], const char mode[])
{
  return fopen(filename, mode);
}

int make_directory(const char path[])
{
  char* buf= (char*) malloc(sizeof(char)*(strlen(path)+1));

  char* p= buf;
  char const * c= path;
  struct stat st;

  do {
    if(*c == '/') {
      *p= '\0';
    
      if(stat(buf, &st)) {
	if(!mkdir(buf, S_IRWXU))
	  cout << "New directory created: " << buf << "\n";
	else {
	  cerr << "Error: Failed to make new directory: " << buf << endl;
	  return false;
	}
      }
    }
  } while((*(p++) = *(c++)));

  free(buf);
  return true;
}

/*
FILE* open(const char filename[], const char mode[])
{
  char* buf= (char*) malloc(sizeof(char)*(strlen(filename)+1));
  //const int n= strlen(filename)+1;  for(int i=0; i<n; ++i) buf[i]= 0; //debug

  char* p= buf;
  char const * c= filename;
  struct stat st;

  do {
    if(*c == '/') {
      *p= '\0';
    
      if(stat(buf, &st)) {
	//cout << "Directory: " << buf << " does not exist\n";
	if(!mkdir(buf, S_IRWXU))
	  cout << "New directory created: " << buf << "\n";
	else {
	  cerr << "Error: Failed to make new directory: " << buf << endl;
	  return 0;
	}
      }
    }
  } while((*(p++) = *(c++)));

  free(buf);
  return fopen(filename, mode);
}  

*/
