program pi_f_mpi
  use  mpi
  ! define variables
  integer          :: i,n,pn,ierr,rank,size
  double precision :: s,s_root,pi
  n = 1000000000
  s = 0.0d0
  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  ! split iterations on processes
  pn = n/size;
  ! compute pi
  do i = 0, pn
    j = i+pn*rank
    s = s + 1.0d0/((4.0d0*j+1.0d0)*(4.0d0*j+3.0d0))
  enddo
  ! reduce sum to root process
  call MPI_REDUCE(s,s_root,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr);
  ! print results (only root process)
  if(rank==0) then
    pi = 8.0d0*s_root
    print *, "approximation to pi =",pi
  endif
  ! end processes
  call MPI_FINALIZE(ierr)
end program pi_f_mpi
