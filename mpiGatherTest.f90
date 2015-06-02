program mpiGatherTest
    ! import mpi
    use mpi
    ! declare variables
    implicit none
    integer              :: i,ierr,nprc,myid,nrow,ncol
    integer, allocatable :: arr(:,:)
    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)
    ! allocate array
    nrow = 4
    ncol = 4
    if (myid == 0) then
        allocate(arr(ncol,nrow*nprc))
    else
        allocate(arr(ncol,nrow))
    endif
    ! fill arrays
    arr(:,:) = myid
    ! gather all values to root array
    call MPI_GATHER(arr,nrow*ncol,MPI_INTEGER,arr,nrow*ncol,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    ! print output
    if (myid == 0) then
        print*,'gathered root array'
        do i = 1, nrow*nprc
            write(*,'(4i4)') arr(:,i)
        enddo
    endif
    ! deallocate array
    deallocate(arr)
    ! end processes
    call MPI_FINALIZE(ierr)
end program mpiGatherTest