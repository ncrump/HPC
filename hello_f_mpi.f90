program hello_f_mpi
  use mpi
  ! define variables
  implicit none
  integer :: ierr,rank,size
  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
  ! print message
  print*,"hello from process",rank
  ! end processes
  call MPI_FINALIZE(ierr)
end program hello_f_mpi
