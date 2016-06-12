#include <stdio.h>

int main(int argc,char **argv){
  // define variables
  int myid,omp_get_thread_num();
  // print threads
  #pragma omp parallel private(myid)
  {
  myid = omp_get_thread_num();
  printf("hello from thread %d\n",myid);
  }
}
