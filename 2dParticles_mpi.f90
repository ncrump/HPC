! Created on Sat Apr 04 18:21:52 2015
! CSI 702, Assignment 4
! Nick Crump

! ** MPI version **

! Subdivides 2D domain along y-axis into blocks for each process
! Assumes 4 processes with 4x1 domain decomposition
! Particle positions and velocities calculated in each block
! Outputs values to text file
program particles
    ! import mpi
    use mpi
    ! declare variables
    implicit none
    integer              :: xymin,xymax,nprt,rprt
    integer              :: i,ierr,nprc,myid
    real(8)              :: dy,dt,cfl,vmax
    real(8), allocatable :: pxn(:),pyn(:),vxn(:),vyn(:)
    real(8), allocatable :: pxr(:),pyr(:),vxr(:),vyr(:)
    ! setup output files
    open(10,file='particles.txt')
    ! set input parameters
    xymin = -4.0
    xymax =  4.0
    nprt  =  100
    cfl   =  1.0
    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)
    ! check domain partition (assumes nprc=4, grid=4x1)
    if (nprc /= 4 .or. mod(xymax-xymin,nprc) /= 0) then
        if (myid == 0) then
            ! if not print error message
            print*, 'warning: sub-domain not divided evenly, check input parameters'
            print*, 'exiting program ...'
            endif
        ! exit
        call MPI_FINALIZE(ierr)
        stop
    endif
    ! allocate arrays
    rprt = int(nprt*nprc)
    allocate(pxn(nprt),pyn(nprt),vxn(nprt),vyn(nprt))
    allocate(pxr(rprt),pyr(rprt),vxr(rprt),vyr(rprt))
    ! partition domain
    dy = (xymax-xymin)/nprc
    ! initialize particle positions
    call randnum(pxn,pyn,nprt,myid)
    ! shift xy-positions to sub-domains
    pxn(:) = (xymax-xymin)*pxn(:) + xymin
    pyn(:) = dy*pyn(:) + xymin + dy*myid
    ! initialize particle velocities
    call grad2d(pxn,pyn,vxn,vyn,nprt)
    ! gather all values to root array
    call MPI_GATHER(pxn,nprt,MPI_DOUBLE_PRECISION,pxr,nprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(pyn,nprt,MPI_DOUBLE_PRECISION,pyr,nprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(vxn,nprt,MPI_DOUBLE_PRECISION,vxr,nprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(vyn,nprt,MPI_DOUBLE_PRECISION,vyr,nprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    ! output values
    if (myid == 0) then
        ! compute allowable time step
        vmax = maxval((vxr(:)**2 + vyr(:)**2)**0.5)
        dt = cfl*(dy/vmax)
        ! print to screen
        print*, 'dt         =',dt
        print*, 'processors =',nprc
        print*, 'particles  =',rprt
        print*, 'xy-domain  =',xymin,xymax
        print*, 'dy         =',dy
        ! print to file
        do i = 1, rprt
            write(10,*) pxr(i),pyr(i),vxr(i),vyr(i)
        enddo
    endif
    ! deallocate arrays
    deallocate(pxn,pyn,vxn,vyn)
    deallocate(pxr,pyr,vxr,vyr)
    ! end processes
    call MPI_FINALIZE(ierr)
end program particles

! subroutine to generate positions
! generates random numbers over [0,1]
subroutine randnum(pxn,pyn,nprt,myid)
    ! declare variables
    implicit none
    integer :: i,nprt,myid,clk
    real(8) :: a,xm,x0,seed
    real(8) :: pxn(nprt),pyn(nprt)
    ! define constants
    a  = 16807.0
    xm = 2147483647.0
    x0 = 2147483711.0
    ! seed number generator
    call system_clock(count=clk)
    seed = clk*(myid+1)
    ! start modulo generator
    do i = 1, nprt
        seed = mod(a*seed,xm)
        pxn(i) = seed/x0
    enddo
    do i = 1, nprt
        seed = mod(a*seed,xm)
        pyn(i) = seed/x0
    enddo
end subroutine randnum

! subroutine to generate velocities
! generates velocities from gradient function
subroutine grad2d(pxn,pyn,vxn,vyn,nprt)
    ! declare variables
    implicit none
    integer :: i,nprt
    real(8) :: xi,yi
    real(8) :: pxn(nprt),pyn(nprt),vxn(nprt),vyn(nprt)
    ! get velocities from gradient function
    do i = 1, nprt
        xi = pxn(i)
        yi = pyn(i)
        vxn(i) = 1.0 + 1.0/(xi**2 + yi**2) - (2*xi**2)/(xi**2 + yi**2)**2
        vyn(i) = -(2*xi*yi)/(xi**2 + yi**2)**2
    enddo
end subroutine grad2d