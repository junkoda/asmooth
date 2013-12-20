#include "mpi_interface.h"
#include <cmath>

using namespace std;

MpiInterface::MpiInterface()
{
  MPI_Comm_size(MPI_COMM_WORLD, &n_node);  
  nc_node= (int) floor(pow(n_node, 1.0/3.0)+0.5);
  if(nc_node*nc_node*nc_node != n_node)
    throw MpiInterfaceError();

  int dims[]= {nc_node, nc_node, nc_node};
  int periods[]= {true, true, true};

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, true, &comm_cart); 
  MPI_Cart_get(comm_cart, 3, dims, periods, coords);
  MPI_Comm_rank(comm_cart, &my_rank);
}

int MpiInterface::sendrecv_to(void* const send_buf, const int nsend,
			      void* const recv_buf, const int nbuf,
			      const int direction_to) const
{
  int dn[3];
  dn[2]= direction_to / 9 - 1;
  dn[1]= (direction_to - 9*(dn[2]+1))/3 - 1;
  dn[0]= direction_to % 3 - 1;
  assert(9*(dn[2]+1) + 3*(dn[1]+1) + (dn[0]+1) == direction_to);

  const int node_from= node_in_direction(-dn[0], -dn[1], -dn[2]);
  const int node_to= node_in_direction(dn[0], dn[1], dn[2]);

  MPI_Status stat;

  //if(my_rank == 0)
  //  fprintf(stderr, "direction %d\n", direction_to);
  MPI_Barrier(comm_cart);
  
  /*
  //debug
  if(nsend > 0) {
    fprintf(stderr, "sendrecv %d -> %d: %d\n", node_to, node_from, direction_to);
    fprintf(stderr, "direction  %d %d %d\n", dn[0], dn[1], dn[2]);
    fprintf(stderr, "coord      %d %d %d\n", coords[0], coords[1], coords[2]);
    int c[]={coords[0], coords[1], coords[2]};
    int rank;
    MPI_Cart_rank(comm_cart, c, &rank);
    fprintf(stderr, "test1 %d\n", rank);
  }
  */
  // C++ function throws error. No need to check return values, i guess...

  int nsend1= nsend, nrecv1;
  int ret1=
    MPI_Sendrecv(&nsend1, 1, MPI_INT, node_to, direction_to /* tag */,
	         &nrecv1, 1, MPI_INT, node_from, direction_to,
	         comm_cart, &stat);
  assert(ret1 == MPI_SUCCESS);

  // debug ***
  //fprintf(stderr, "%d: %d %d %d %d < %d\n", my_rank, node_from, node_to, nsend, nrecv1, nbuf);
  //MPI_Barrier(comm_cart);

  int ret=
    MPI_Sendrecv(send_buf, nsend,  MPI_FLOAT, node_to, direction_to /* tag */,
	         recv_buf, nrecv1, MPI_FLOAT, node_from, direction_to,
	         comm_cart, &stat);
  
  //  MPI_Sendrecv(send_buf, nsend, MPI_FLOAT, node_to, direction_to /* tag */,
  //recv_buf, nbuf, MPI_FLOAT, node_from, direction_to,
  //comm_cart, &stat);

  if(ret == MPI_ERR_TRUNCATE)
    abort("Not enough space to sendrecv buffer particles\n");
  assert(ret == MPI_SUCCESS);

  int nrecv= 0;
  MPI_Get_count(&stat, MPI_FLOAT, &nrecv);

  return nrecv; // return number of recieved particles
}

int MpiInterface::bcast(void *buf, int count, MPI_Datatype datatype) const
{
  return MPI_Bcast(buf, count, datatype, 0, comm_cart);
}
