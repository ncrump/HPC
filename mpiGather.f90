program mpiGather
    ! import mpi
    use mpi
    ! declare variables
    implicit none
    integer              :: i,ierr,nprc,myid,nval,ntot
    integer, allocatable :: myarr(:),rootarr(:)
    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)
    ! allocate arrays
    nval = 2
    ntot = int(nval*nprc)
    allocate(myarr(nval),rootarr(ntot))
    ! loop through array
    print*,'processor array'
    do i = 1, nval
        myarr(i) = i*myid
    enddo
    print*,myid,myarr
    ! gather all values to root array
    call MPI_GATHER(myarr,nval,MPI_INTEGER,  &
                    rootarr,nval,MPI_INTEGER,&
                    0,MPI_COMM_WORLD,ierr)
    ! print output
    if (myid == 0) then
        print*,'gathered root array'
        do i = 1, ntot
            print*,rootarr(i)
        enddo
    endif
    ! deallocate arrays
    deallocate(myarr,rootarr)
    ! end processes
    call MPI_FINALIZE(ierr)
end program mpiGather