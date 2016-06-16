#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


// line model
void line_model(double a,double b,double *xmod,double *ymod,
                double *yerr,int npts){
  for(int i=0; i<npts; i++){
    xmod[i] = i;
    ymod[i] = a*xmod[i]+b;
    yerr[i] = fabs(ymod[i])*0.1;
  }
}

// log posterior
double ln_posterior(double ak,double bk,double *xmod,double *ymod,double *yerr,
                    double beta2,int npts){
  double ln_prior = log(beta2+ak*ak)+log(beta2+bk*bk);
  double ln_likel = 0.0;
  for(int j=0; j<npts; j++){
    double ydff = ymod[j]-(ak*xmod[j]+bk);
    ydff = (ydff*ydff)/(yerr[j]*yerr[j]);
    ln_likel += ydff;
  }
  return -ln_prior-ln_likel;
}

// permute index
void permute(int *rindex,int walkers,int seed){
  srand(time(NULL)*seed);
  for(int i=0; i<walkers-1; i++){
    int j = i+rand()/(RAND_MAX/(walkers-i)+1);
    int tt = rindex[j];
    rindex[j] = rindex[i];
    rindex[i] = tt;
  }
}


int main(int argc,char **argv){

  // inputs
  // ----------------------
  // set size
  int    walkers = 1000;
  int    steps   = 10000;
  double a       = 5.6;  // line slope
  double b       = 2.2;  // line intercept
  int    npts    = 100;  // number of points
  int    npar    = 2;    // number of parameters
  double beta    = 0.2;  // hyperparameter
  // ----------------------

  // initialize MPI
  int rank,size;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  // split walkers into groups for processes
  int pwalkers = walkers/(float)size;

  // allocate arrays
  double *xmod   = malloc(npts*sizeof(double));
  double *ymod   = malloc(npts*sizeof(double));
  double *yerr   = malloc(npts*sizeof(double));
  double *yfit   = malloc(npts*sizeof(double));
  double *post0  = malloc(pwalkers*sizeof(double));
  double *wlkrs1 = malloc(pwalkers*sizeof(double));
  double *wlkrs2 = malloc(pwalkers*sizeof(double));
  double **rran  = malloc(steps*sizeof(double));
  double **zham  = malloc(steps*sizeof(double));
  for(int i=0; i<steps; i++){
    rran[i] = malloc(pwalkers*sizeof(double));
    zham[i] = malloc(pwalkers*sizeof(double));
  }
  int    *rindex  = malloc(walkers*sizeof(int));
  double *rpost0  = malloc(walkers*sizeof(double));
  double *rwlkrs1 = malloc(walkers*sizeof(double));
  double *rwlkrs2 = malloc(walkers*sizeof(double));
  double *rpost0_tmp  = malloc(walkers*sizeof(double));
  double *rwlkrs1_tmp = malloc(walkers*sizeof(double));
  double *rwlkrs2_tmp = malloc(walkers*sizeof(double));

  // get line model
  line_model(a,b,xmod,ymod,yerr,npts);

  // seed random numbers
  srand(time(NULL)*(rank+1));

  // initialize walkers (only root process)
  if(rank==0){
    for(int i=0; i<walkers; i++){
      rwlkrs1[i] = rand()/(double)RAND_MAX;
      rwlkrs2[i] = rand()/(double)RAND_MAX;
      rindex[i] = i;
    }
  }

  // scatter the chunks from root to others
  MPI_Scatter(rwlkrs1,pwalkers,MPI_DOUBLE,wlkrs1,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Scatter(rwlkrs2,pwalkers,MPI_DOUBLE,wlkrs2,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // initialize acceptance
  double accept  = 0.0;
  double raccept = 0.0;

  // initialize mcmc hammer
  double aham = 2.0;
  double zmin = sqrt(1.0/aham);
  double zmax = sqrt(aham);
  double zdff = zmax-zmin;
  int    phalf = pwalkers/2.0;
  // get z values
  for(int i=0; i<steps; i++){
    for(int j=0; j<pwalkers; j++){
      double randm = rand();
      randm = randm/RAND_MAX;
      randm = randm*zdff+zmin;
      randm = randm*randm;
      zham[i][j] = randm;
    }
  }

  // initialize posterior
  double beta2 = beta*beta;
  for(int i=0; i<pwalkers; i++){
    double ak = wlkrs1[i];
    double bk = wlkrs2[i];
    post0[i] = ln_posterior(ak,bk,xmod,ymod,yerr,beta2,npts);
  }

  // gather the chunks from others to root
  MPI_Gather(post0,pwalkers,MPI_DOUBLE,rpost0,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);

  // generate random numbers
  for(int i=0; i<steps; i++){
    for(int j=0; j<pwalkers; j++){
      double randm = rand();
      randm = randm/RAND_MAX;
      rran[i][j] = randm;
    }
  }

  // begin mcmc
  int i,j,k,w1,w2,sign;
  for(i=0; i<steps; i++){

    // loop over sets of half-walkers
    for(j=0; j<2; j++){

      // split walkers into half-groups
      if(j==0){
        // only root permutes index
        if(rank==0){
          permute(rindex,walkers,(i+1)*(j+1));
          // root shuffles walkers
          for(int t=0; t<walkers; t++){
            int tndx = rindex[t];
            rwlkrs1_tmp[t] = rwlkrs1[tndx];
            rwlkrs2_tmp[t] = rwlkrs2[tndx];
            rpost0_tmp[t]  = rpost0[tndx];
          }
        }
        // scatter the chunks from root to others
        MPI_Scatter(rwlkrs1_tmp,pwalkers,MPI_DOUBLE,wlkrs1,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Scatter(rwlkrs2_tmp,pwalkers,MPI_DOUBLE,wlkrs2,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Scatter(rpost0_tmp,pwalkers,MPI_DOUBLE,post0,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
        // set group 1
        w1 = 0;
        w2 = phalf;
        sign = 1;
      }
      if(j==1){
        // set group 2
        w1 = phalf;
        w2 = pwalkers;
        sign = -1;
      }

      // loop over half-walkers
      for(k=w1; k<w2; k++){

        // get walker k from group 1
        int wk = k;
        double ak0 = wlkrs1[wk];
        double bk0 = wlkrs2[wk];

        // get walker j from group 2
        int wj = k+sign*phalf;
        double aj0 = wlkrs1[wj];
        double bj0 = wlkrs2[wj];

        // do stretch move
        double ak1 = aj0 + zham[i][wk]*(ak0-aj0);
        double bk1 = bj0 + zham[i][wk]*(bk0-bj0);

        // get new posterior
        double post1 = ln_posterior(ak1,bk1,xmod,ymod,yerr,beta2,npts);

        // accept/reject move
        double prob = post1-post0[wk]+(npts-1)*log(zham[i][wk]);
        // accept
        if(prob>0){
          wlkrs1[wk] = ak1;
          wlkrs2[wk] = bk1;
          post0[wk] = post1;
          accept += 1;
        }
        // otherwise accept with some probability
        else{
          if(prob>log(rran[i][wk])){
            wlkrs1[wk] = ak1;
            wlkrs2[wk] = bk1;
            post0[wk] = post1;
            accept += 1;
          }
          // otherwise reject
        }
      }
    }

    // gather the chunks from others to root
    MPI_Gather(wlkrs1,pwalkers,MPI_DOUBLE,rwlkrs1,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(wlkrs2,pwalkers,MPI_DOUBLE,rwlkrs2,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Gather(post0,pwalkers,MPI_DOUBLE,rpost0,pwalkers,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }

  // reduce acceptance to root process
  MPI_Reduce(&accept,&raccept,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  // get final values (only root process)
  if(rank==0){
    // get average of walkers
    double avg_a = 0.0;
    double avg_b = 0.0;
    for(int i=0; i<walkers; i++){
      avg_a += rwlkrs1[i];
      avg_b += rwlkrs2[i];
    }
    avg_a = avg_a/walkers;
    avg_b = avg_b/walkers;

    // get final fit
    double chisq = 0.0;
    for(int i=0; i<npts; i++){
      yfit[i] = avg_a*xmod[i]+avg_b;
      double chi = (ymod[i]-yfit[i])/yerr[i];
      chisq += chi*chi;
    }
    chisq = chisq/(npts-npar);
    raccept = raccept/(walkers*steps);

    // print results
    printf("walkers      = %i\n",walkers);
    printf("steps        = %i\n",steps);
    printf("acceptance   = %f\n",raccept);
    printf("chisq        = %f\n",chisq);
    printf("model params = %f %f\n",a,b);
    printf("fit params   = %f %f\n",avg_a,avg_b);
    printf("\n");
  }

  // free arrays
  free(xmod);
  free(ymod);
  free(yerr);
  free(yfit);
  free(post0);
  free(wlkrs1);
  free(wlkrs2);
  free(rran);
  free(zham);
  free(rindex);
  free(rpost0);
  free(rwlkrs1);
  free(rwlkrs2);
  free(rpost0_tmp);
  free(rwlkrs1_tmp);
  free(rwlkrs2_tmp);

  // end processes
  MPI_Finalize();

}
