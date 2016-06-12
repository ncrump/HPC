#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc,char **argv){

  // define variables
  int    i,j;
  double msd,xij,yij,randm;

  // inputs
  // ----------------------
  // set size
  int walkers = 2000;
  int steps   = 50000;
  // ----------------------

  // allocate arrays
  double *xpos = malloc(walkers*sizeof(double));
  double *ypos = malloc(walkers*sizeof(double));
  double *thet = malloc(walkers*steps*sizeof(double));
  int    *step = malloc(steps*sizeof(int));
  double *rmsd = malloc(steps*sizeof(double));

  // initialize
  double dr = 0.2;
  double pi = 4.0*atan(1.0);
  for(i=0; i<walkers; i++){
    xpos[i] = 0.0;
    ypos[i] = 0.0;
  }

  // generate random angles
  srand(time(NULL));
  for(i=0; i<steps; i++){
    for(j=0; j<walkers; j++){
      randm = rand();
      randm = (randm/RAND_MAX)*2.0*pi;
      thet[i*walkers+j] = randm;
    }
  }

  // random walk
  #pragma omp parallel private(i,j,xij,yij)
  for(i=0; i<steps; i++){
    msd = 0.0;
    #pragma omp barrier
    #pragma omp for reduction(+:msd)
    for(j=0; j<walkers; j++){
      xpos[j] += dr*cos(thet[i*walkers+j]);
      ypos[j] += dr*sin(thet[i*walkers+j]);
      xij = xpos[j];
      yij = ypos[j];
      // get displacement
      msd += xij*xij + yij*yij;
    }
    // store values to array
    #pragma omp single
    step[i] = i+1;
    #pragma omp single
    rmsd[i] = sqrt(msd/walkers);
  }

  // write output to file
  FILE *f = fopen("random_walk_c_omp.txt","w");
  for(i=0; i<steps; i++){
    fprintf(f,"%i  %f\n",step[i],rmsd[i]);
  }
  fclose(f);

  // free arrays
  free(xpos);
  free(ypos);
  free(thet);
  free(step);
  free(rmsd);

}
