#include <stdio.h>
#include <mpi.h>

int main(int argc,char **argv){
  // define variables
  int rank,size;
  // initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  // print message
  printf("hello from process %d\n",rank);
  // end processes
  MPI_Finalize();
}
