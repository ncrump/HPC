#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int main(int argc,char **argv){

  // define variables
  int    i,j,pwalkers,rank,size;
  double msd,msd_root,xij,yij,randm;

  // inputs
  // ----------------------
  // set size
  int walkers = 2000;
  int steps   = 50000;
  // ----------------------

  // initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // split walkers into groups for processes
  pwalkers = walkers/(float)size;

  // allocate arrays
  double *xpos = malloc(pwalkers*sizeof(double));
  double *ypos = malloc(pwalkers*sizeof(double));
  double *thet = malloc(pwalkers*steps*sizeof(double));
  int    *step = malloc(steps*sizeof(int));
  double *rmsd = malloc(steps*sizeof(double));

  // initialize
  double dr = 0.2;
  double pi = 4.0*atan(1.0);
  for(i=0; i<pwalkers; i++){
    xpos[i] = 0.0;
    ypos[i] = 0.0;
  }

  // generate random angles
  srand(time(NULL)*rank+1);
  for(i=0; i<steps; i++){
    for(j=0; j<pwalkers; j++){
      randm = rand();
      randm = (randm/RAND_MAX)*2.0*pi;
      thet[i*pwalkers+j] = randm;
    }
  }

  // random walk

  for(i=0; i<steps; i++){
    msd = 0.0;
    for(j=0; j<pwalkers; j++){
      xpos[j] += dr*cos(thet[i*pwalkers+j]);
      ypos[j] += dr*sin(thet[i*pwalkers+j]);
      xij = xpos[j];
      yij = ypos[j];
      // get displacement
      msd += xij*xij + yij*yij;
    }
    // reduce msd to root process
    MPI_Reduce(&msd,&msd_root,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    // store values to arrays (only root process)
    if(rank==0){
      step[i] = i+1;
      rmsd[i] = sqrt(msd_root/pwalkers);
    }
  }

  // output to file (only root process)
  if(rank==0){
    FILE *f = fopen("random_walk_c_mpi.txt","w");
    for(i=0; i<steps; i++){
      fprintf(f,"%i %f\n",step[i],rmsd[i]);
    }
    fclose(f);
  }

  // free arrays
  free(xpos);
  free(ypos);
  free(thet);
  free(step);
  free(rmsd);

  // end processes
  MPI_Finalize();

}
