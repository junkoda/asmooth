#ifndef K_NEIGHBORS_H
#define K_NEIGHBORS_H 1

#include <limits>

//
// CLASS KthValue
// 
class KNeighbors {
 public:
  KNeighbors() {
    const float max_float= std::numeric_limits<float>::max();
    x= new float[k+1];
    s= new Particle*[k+1];
    for(int i=0; i<k; ++i) {
      x[i]= max_float;
      s[i]= 0;
    }
    x[k]= -max_float;
  }
  ~KNeighbors() {
    delete [] x;
    delete [] s;
  }
  void push(const float val, Particle * const p) {
    if(val < x[0]) {
      int i;
      for(i=0; val < x[i+1]; ++i) {
	x[i]= x[i+1];
	s[i]= s[i+1];
      }
      x[i]= val;
      s[i]= p;
    }
  }
  float get_kth_value() const {
    return x[0];
  }
  void clear() {
    for(int i=0; i<k; ++i) {
      x[i]= std::numeric_limits<float>::max();
      s[i]= 0;
    }
  }
  Particle* nbr(const int i){ return s[i]; }
  float r2(const int i){ return x[i]; }
  static const int k= 33;
 private:
  float* x;
  Particle** s;
};

#endif
