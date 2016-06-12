#include <stdio.h>

int main(int argc,char **argv){
  // define variables
  int    i,n = 1000000000;
  double pi,s = 0.0;
  // compute pi
  for(i=0; i<n; i++){
    s += 1.0/((4.0*i+1.0)*(4.0*i+3.0));
  }
  pi = 8.0*s;
  // print results
  printf("approximation to pi = %1.15f\n",pi);
}
