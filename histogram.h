#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <assert.h>

class Histogram {
 public:
  Histogram();
  Histogram(const float xmin, const float xmax, const int nbin);
  ~Histogram();
  void initialize(const float xmin, const float xmax, const int nbin);
  void clear();
  long long ndata() const { return ntot; }
  int operator[](const int i) const { return n[i]; }
  void add(const float x) {

#ifdef _OPENMP
#pragma omp atomic
#endif
    ntot++;

    if(xmin <= x && x < xmax) {
      int i= (int)((x-xmin)/(xmax-xmin)*nbin);
      assert(0<=i && i<nbin);

#ifdef _OPENMP
#pragma omp atomic
#endif
      n[i]++;

    }
  }
  void write(FILE* const fp) const;
  void read(FILE* const fp);
  void write_header(FILE* const fp) const;
  void print() const;
  Histogram& operator+=(const Histogram&);
 private:
  int *n;
  long long ntot;
  double sum;
  float xmin, xmax;
  int nbin;
};


#endif
