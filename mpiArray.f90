program mpiArray
    ! import mpi
    use mpi
    ! declare variables
    implicit none
    integer              :: i,ierr,nprc,myid
    real(8), allocatable :: x(:)
    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)
    ! allocate array
    allocate(x(nprc))
    ! loop through array
    do i = 1, nprc
        x(i) = myid
    enddo
    ! print output
    print*, myid,x
    ! deallocate array
    deallocate(x)
    ! end processes
    call MPI_FINALIZE(ierr)
end program mpiArray