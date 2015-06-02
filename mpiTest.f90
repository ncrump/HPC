program mpiTest
    use mpi
    implicit none
    integer :: ierr,nprc,myid,i

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)
    do i=1,2
        print*, 'proc=',myid, 'iter=',i
    enddo
    call MPI_FINALIZE(ierr)
end program mpiTest
