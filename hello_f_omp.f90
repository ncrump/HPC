program hello_f_omp
  ! define variables
  implicit none
  integer :: myid,omp_get_thread_num
  ! print threads
  !$omp parallel private(myid)
  myid = omp_get_thread_num()
  print*,"hello from thread",myid
  !$omp end parallel
end program hello_f_omp
