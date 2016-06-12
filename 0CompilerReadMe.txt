Compile Options
--------------------------------
****Standard****
gfortran mycode.f90 -o mycode
gcc      mycode.c   -o mycode

****OpenMP****
gfortran mycode.f90 -fopenmp -o mycode
gcc      mycode.f90 -fopenmp -o mycode

****MPI****
mpif90 mycode.f90 -o mycode
mpicc  mycode.c   -o mycode
--------------------------------

Run Options
-------------------------------
****Standard/OpenMP****
./mycode

****Timed****
time ./mycode

****MPI****
mpirun -np 4 ./mycode
-------------------------------