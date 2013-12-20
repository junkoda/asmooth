#include "histogram.h"
#include <iostream>

using namespace std;

Histogram::Histogram() :
  n(0)
{

}

Histogram::Histogram(const float xmin_, const float xmax_, const int nbin_) :
  xmin(xmin_), xmax(xmax_), nbin(nbin_)
{
  n= (int*) malloc(sizeof(int)*nbin);
  assert(n);
}

void Histogram::initialize(const float xmin_, const float xmax_, const int nbin_)
{
  xmin= xmin_; xmax= xmax_; nbin= nbin_;
  n= (int*) malloc(sizeof(int)*nbin);
  assert(n);
}

Histogram::~Histogram()
{
  free(n);
}

void Histogram::clear()
{
  ntot=0;
  for(int j=0; j<nbin; ++j)
    n[j]= 0;
}

void Histogram::write_header(FILE* const fp) const
{
  fwrite(&xmin, sizeof(float), 1, fp);
  fwrite(&xmax, sizeof(float), 1, fp);
  fwrite(&nbin, sizeof(int), 1, fp);

}

void Histogram::write(FILE* const fp) const
{
  float fsum= (float) sum;
  //int ntot_int= (int) ntot;
  fwrite(&ntot, sizeof(long long), 1, fp);
  assert(sizeof(long long) == 8);
  //cerr << "ntot " << ntot << "\n";
  fwrite(&fsum, sizeof(float), 1, fp);
  fwrite(n, sizeof(int), nbin, fp);
}

void Histogram::print() const
{
  printf("range %e %e %d\n", xmin, xmax, nbin);

  const float dx= (xmax - xmin)/nbin;
  for(int j=0; j<nbin; ++j) {
    float xmid= xmin + (j+0.5f)*(xmax - xmin)/nbin;
    if(n[j] > 0)
      printf("%e %e\n", xmid, double(n[j])/ntot/dx);
  }
}

Histogram& Histogram::operator+=(const Histogram& h)
{
  assert(n);
  assert(nbin == h.nbin);
  assert(xmin == h.xmin);
  assert(xmax == h.xmax);

  ntot += h.ntot;
  sum += h.sum;

  for(int j=0; j<nbin; ++j)
    n[j] += h.n[j];

  return *this;
}

void Histogram::read(FILE* const fp)
{
  int ntot_int= 0;
  assert(n);
  //fread(&ntot_int, sizeof(int), 1, fp);
  //ntot= ntot_int; //debug (old)
  fread(&ntot, sizeof(long long), 1, fp);
  fread(&sum, sizeof(float), 1, fp);
  fread(n, sizeof(int), nbin, fp);
}
