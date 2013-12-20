#ifndef MPI_INTERFACE_H
#define MPI_INTERFACE_H 1
#include <mpi.h>
#include <iostream>
#include <ostream>
// #include <cmath>
#include <assert.h>

//
// Interface for Message Passing Interface
//

class MpiInterface {
 public:
  MpiInterface();
  int rank() const { return my_rank; }
  int const * coordinate() const {
    return coords;
  }
  int index() const {
    return coords[0]
         + coords[1]*nc_node
         + coords[2]*nc_node*nc_node;
  }
  int size() const {
    return n_node;
  }
  int nc() const {
    return nc_node;
  }
  void abort(const char msg[]) const {
    std::cerr << msg;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  void print(const char msg[]) const {
    if(index() == 0)
       std::cout << msg;
  }
  int sendrecv_to(void * const buf_send, const int nsend,
		  void* const buf_recv, const int nrecv, 
		   const int direction_to) const;
  int node_in_direction(const int dx, const int dy, const int dz) const {
    int c[]= {(coords[0]+dx+nc_node) % nc_node, 
	      (coords[1]+dy+nc_node) % nc_node,
	      (coords[2]+dz+nc_node) % nc_node};
    int rank= -1;
    MPI_Cart_rank(comm_cart, c, &rank);
    return rank;
  }

  int bcast(void *buf, int count, MPI_Datatype datatype) const;

  struct MpiInterfaceError{};
 private:
  int my_index, n_node, nc_node, my_rank;
  int coords[3];
  MPI_Comm comm_cart;
};

#endif

