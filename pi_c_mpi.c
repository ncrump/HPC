#include <stdio.h>
#include <mpi.h>

int main(int argc,char **argv){
  // define variables
  int    i,j,rank,size,pn,n = 1000000000;
  double pi,s_root,s = 0.0;
  // initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  // split iterations on processes
  pn = n/(float)size;
  // compute pi
  for(i=0; i<pn; i++){
    j = i+pn*rank;
    s += 1.0/((4.0*j+1.0)*(4.0*j+3.0));
  }
  // reduce sum to root process
  MPI_Reduce(&s,&s_root,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  // print results (only root process)
  if(rank==0){
    pi = 8.0*s_root;
    printf("approximation to pi = %1.15f\n",pi);
  }
  // end processes
  MPI_Finalize();
}
