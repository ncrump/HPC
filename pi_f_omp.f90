program pi_f_omp
  ! define variables
  integer          :: i,n
  double precision :: s,pi
  n = 1000000000
  s = 0.0d0
  ! compute pi
  !$omp parallel do private(i) shared(n) reduction(+:s)
  do i = 0, n
    s = s + 1.0d0/((4.0d0*i+1.0d0)*(4.0d0*i+3.0d0))
  enddo
  !$omp end parallel do
  pi = 8.0d0*s
  ! print results
  print *, "approximation to pi =",pi
end program pi_f_omp
