#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// line model
void line_model(double a,double b,double *xmod,double *ymod,
                double *yerr,int npts){
  for(int i=0; i<npts; i++){
    xmod[i] = i;
    ymod[i] = a*xmod[i]+b;
    yerr[i] = abs(ymod[i])*0.1;
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

void permute(int *index,int walkers,int seed){
  srand(time(NULL)*seed);
  for(int i=0; i<walkers-1; i++){
    int j = i+rand()/(RAND_MAX/(walkers-i)+1);
    int t = index[j];
    index[j] = index[i];
    index[i] = t;
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

  // allocate arrays
  double *xmod   = malloc(npts*sizeof(double));
  double *ymod   = malloc(npts*sizeof(double));
  double *yerr   = malloc(npts*sizeof(double));
  double *yfit   = malloc(npts*sizeof(double));
  int    *index  = malloc(walkers*sizeof(int));
  double *post0  = malloc(walkers*sizeof(double));
  double **wlkrs = malloc(walkers*sizeof(double));
  for(int i=0; i<walkers; i++){
    wlkrs[i] = malloc(npar*sizeof(double));
  }
  double **rran = malloc(steps*sizeof(double));
  double **zham = malloc(steps*sizeof(double));
  for(int i=0; i<steps; i++){
    rran[i] = malloc(walkers*sizeof(double));
    zham[i] = malloc(walkers*sizeof(double));
  }

  // get line model
  line_model(a,b,xmod,ymod,yerr,npts);

  // initialize walkers
  srand(time(NULL));
  for(int i=0; i<walkers; i++){
    wlkrs[i][0] = rand()/(double)RAND_MAX;
    wlkrs[i][1] = rand()/(double)RAND_MAX;
    index[i]    = i;
  }

  // initialize acceptance
  double accept = 0.0;

  // initialize mcmc hammer
  double aham = 2.0;
  double zmin = sqrt(1.0/aham);
  double zmax = sqrt(aham);
  double zdff = zmax-zmin;
  int    half = walkers/2.0;
  // get z values
  for(int i=0; i<steps; i++){
    for(int j=0; j<walkers; j++){
      double randm = rand();
      randm = randm/RAND_MAX;
      randm = randm*zdff+zmin;
      randm = randm*randm;
      zham[i][j] = randm;
    }
  }

  // initialize posterior
  double beta2 = beta*beta;
  for(int i=0; i<walkers; i++){
    double ak = wlkrs[i][0];
    double bk = wlkrs[i][1];
    post0[i] = ln_posterior(ak,bk,xmod,ymod,yerr,beta2,npts);
  }

  // generate random numbers
  for(int i=0; i<steps; i++){
    for(int j=0; j<walkers; j++){
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
        // permute index
        permute(index,walkers,(i+1)*(j+1));
        // set group 1
        w1 = 0;
        w2 = half;
        sign = 1;
      }
      if(j==1){
        // set group 2
        w1 = half;
        w2 = walkers;
        sign = -1;
      }

      // loop over half-walkers
      #pragma omp parallel for private(k)
      for(k=w1; k<w2; k++){

        // get walker k from group 1
        int wk = index[k];
        double ak0 = wlkrs[wk][0];
        double bk0 = wlkrs[wk][1];

        // get walker j from group 2
        int wj = index[k+sign*half];
        double aj0 = wlkrs[wj][0];
        double bj0 = wlkrs[wj][1];

        // do stretch move
        double ak1 = aj0 + zham[i][wk]*(ak0-aj0);
        double bk1 = bj0 + zham[i][wk]*(bk0-bj0);

        // get new posterior
        double post1 = ln_posterior(ak1,bk1,xmod,ymod,yerr,beta2,npts);

        // accept/reject move
        double prob = post1-post0[wk]+(npts-1)*log(zham[i][wk]);
        // accept
        if(prob>0){
          wlkrs[wk][0] = ak1;
          wlkrs[wk][1] = bk1;
          post0[wk] = post1;
          accept += 1;
        }
        // otherwise accept with some probability
        else{
          if(prob>log(rran[i][wk])){
            wlkrs[wk][0] = ak1;
            wlkrs[wk][1] = bk1;
            post0[wk] = post1;
            accept += 1;
          }
          // otherwise reject
        }
      }
    }
  }

  // get average of walkers
  double avg_a = 0.0;
  double avg_b = 0.0;
  for(int i=0; i<walkers; i++){
    avg_a += wlkrs[i][0];
    avg_b += wlkrs[i][1];
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
  accept = accept/(walkers*steps);

  // print results
  printf("walkers      = %i\n",walkers);
  printf("steps        = %i\n",steps);
  printf("acceptance   = %f\n",accept);
  printf("chisq        = %f\n",chisq);
  printf("model params = %f %f\n",a,b);
  printf("fit params   = %f %f\n",avg_a,avg_b);
  printf("\n");

  // free arrays
  free(xmod);
  free(ymod);
  free(yerr);
  free(yfit);
  free(index);
  free(post0);
  free(wlkrs);
  free(rran);
  free(zham);

}
