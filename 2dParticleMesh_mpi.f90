! Created on Sat May 09 14:07:31 2015
! CSI 702, Assignment 6
! Nick Crump

! ** MPI version **  ** warning: this is a hack in progress...**

! 2D parallel particle mesh simulation for electric charges
! Subdivides domain along y-axis into blocks for each process
! Initializes particles to random positions and zero velocities
! Interpolates particle charge density to nearest grid points
! Solves Poisson equation using Jacobi/SOR iteration method
! Computes gradient of potential on grid to get E-field
! Interpolates forces on grid to nearest particles
! Advances particles using Verlet method
! Outputs values to text files

program particlemesh

    ! import mpi
    use mpi
    ! declare variables
    implicit none
    character            :: frmt*20,npts*7
    integer              :: ierr,nprc,myid,tag,myprt,mystat(MPI_STATUS_SIZE)
    integer              :: i,j,nx,ny,nysub,iter,nprt,nstp,step,nsav,ndx,ndy,t1,t2,rate
    real(8)              :: pi,xmin,xmax,ymin,ymax,tol,tstp,dt
    real(8)              :: dx,dy,py,relax,myerr,mxerr,Eold,Enew,Lx,Ly
    real(8), allocatable :: xgrid(:),ygrid(:),chrg(:,:),Epot(:,:),Ex(:,:),Ey(:,:)
    real(8), allocatable :: pxn(:),pyn(:),vxn(:),vyn(:),Fxn(:),Fyn(:),qn(:)
    real(8), allocatable :: tmpchrg(:,:),allEx(:,:),allEy(:,:),upper(:),lower(:)
    real(8), allocatable :: pxr(:),pyr(:),vxr(:),vyr(:)

    ! set input parameters
    ! -------------------------
    myprt = 50                 ! number of particles per process
    nstp  = 1000               ! number of time steps
    nsav  = 10                 ! steps between saves
    tol   = 1e-5               ! accuracy tolerance
    xmin  = 0    ; xmax = 20   ! x domain
    ymin  = 0    ; ymax = 20   ! y domain
    nx    = 100  ;  ny  = 100  ! points in domain
    dt    = 0.005              ! time step
    ! -------------------------

    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)

    ! set parameters from input
    pi   = 4.0*atan(1.0)
    tstp = 0
    ! grid domain
    Lx = xmax-xmin
    Ly = ymax-ymin
    ! total number of particles
    nprt = myprt*nprc
    ! set relaxation parameter
    relax = 0.8

    ! get domain partition (y-axis subdivided)
    nysub = ny/nprc
    dx    = Lx/(nx-1)
    dy    = Ly/(ny-1)
    py    = Ly/nprc

    ! check domain is square and processes evenly subdivide y-domain
    if (Lx /= Ly .or. nx /= ny .or. mod(ny,nprc) /= 0) then
        if (myid == 0) then
            ! if not print error message
            print*, 'warning: sub-domain not divided evenly, check input parameters'
            print*, 'exiting program ...'
        endif
        ! exit
        call MPI_FINALIZE(ierr)
        stop
    endif

    ! only root process does this
    if (myid == 0) then
        ! setup output files
        open(10,file='xygrid.txt')
        open(11,file='potential.txt')
        open(12,file='Exfield.txt')
        open(13,file='Eyfield.txt')
        open(14,file='trajectory.txt')
        ! set format write string for output files
        write(npts,'(i7)') nx
        frmt = '('//trim(npts)//'f12.6)'
        ! start timer
        call system_clock(t1,rate)
        ! allocate full arrays to root process
        allocate(xgrid(nx),ygrid(ny),chrg(nx,ny),Epot(nx,ny),Ex(nx,ny),Ey(nx,ny))
        allocate(pxn(myprt),pyn(myprt),vxn(myprt),vyn(myprt),Fxn(myprt),Fyn(myprt),qn(myprt))
        allocate(tmpchrg(nx,ny),allEx(nx,ny),allEy(nx,ny),upper(nx),lower(nx))
        allocate(pxr(nprt),pyr(nprt),vxr(nprt),vyr(nprt))
    ! all other processes do this
    else
       ! allocate sub arrays to all others
       allocate(xgrid(nx),ygrid(nysub),chrg(nx,ny),Epot(nx,nysub),Ex(nx,nysub),Ey(nx,nysub))
       allocate(pxn(myprt),pyn(myprt),vxn(myprt),vyn(myprt),Fxn(myprt),Fyn(myprt),qn(myprt))
       allocate(tmpchrg(nx,ny),allEx(nx,ny),allEy(nx,ny),upper(nx),lower(nx))
    endif

    ! get x-grid
    do i = 1, nx
        xgrid(i) = xmin + dx*(i-1)
    enddo
    ! get y-grid (subdivided)
    do j = 1, nysub
        ygrid(j) = ymin + dy*(j-1) + dy*myid*nysub
    enddo

    ! initialize particle positions and charges
    call randnum(pxn,pyn,qn,myprt,myid)
    pxn(:) = Lx*pxn(:) + xmin
    pyn(:) = py*pyn(:) + ymin + py*myid
    qn(:)  = (2*qn(:)-1)/(dx*dy)

    ! initialize velocities to zero
    vxn(:) = 0.0
    vyn(:) = 0.0

    ! iterate over time steps
    do step = 1, nstp
        tstp = tstp + dt
        iter  = 0
        mxerr = 1

        ! zero arrays
        chrg(:,:)    = 0.0
        Epot(:,:)    = 0.0
        Ex(:,:)      = 0.0
        Ey(:,:)      = 0.0
        Fxn(:)       = 0.0
        Fyn(:)       = 0.0
        tmpchrg(:,:) = 0.0
        allEx(:,:)   = 0.0
        allEy(:,:)   = 0.0

        ! compute charge density on grid from charges on particles
        ! **Note: all processes see the charge grid since particles migrate**
        do i = 1, myprt
            ndx = int(pxn(i)/dx)+1
            ndy = int(pyn(i)/dy)+1
            if (ndx < 1)  ndx = 1
            if (ndy < 1)  ndy = 1
            if (ndx > nx) ndx = nx
            if (ndy > ny) ndy = ny
            tmpchrg(ndx,ndy) = tmpchrg(ndx,ndy) + qn(i)
        enddo

        ! sum all charge grids across processes
        call MPI_AllREDUCE(tmpchrg(:,:),chrg(:,:),nx*ny,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! iterate over potential until error is within tolerance
        do while (tol < mxerr)
            mxerr = 0
            iter  = iter + 1
            do j = 1, nysub
                do i = 2, nx-1
                    ! get bottom boundary value
                    if (j == 1 .and. myid == 0) then
                        Epot(i,j) = 0.0
                    ! get top boundary value
                    elseif (j == nysub .and. myid == nprc-1) then
                        Epot(i,j) = 0.0
                    else
                        Eold = Epot(i,j)
                        ! get interior points of electric potential using SOR method
                        if (j == 1) then
                            Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+lower(i)+chrg(i,j+nysub*myid)*dx*dy)
                        elseif (j == nysub) then
                            Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+upper(i)+Epot(i,j-1)+chrg(i,j+nysub*myid)*dx*dy)
                        else
                            Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+Epot(i,j-1)+chrg(i,j+nysub*myid)*dx*dy)
                        endif
                        ! update potential
                        Epot(i,j) = Enew + relax*(Enew-Eold)
                        ! accumulate instantaneous error
                        myerr = myerr + abs(Enew-Eold)
                    endif
                enddo
            enddo
            ! get max average error from all processes
            myerr = myerr/(nx*nysub)
            call MPI_AllREDUCE(myerr,mxerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
            ! get process communication
            tag = iter
            ! send upper subdomain data and receive lower subdomain data
            if (myid < nprc-1) then
                call MPI_SEND(Epot(:,nysub),nx,MPI_DOUBLE_PRECISION,myid+1,tag,MPI_COMM_WORLD,ierr)
                call MPI_RECV(upper(:),nx,MPI_DOUBLE_PRECISION,myid+1,tag,MPI_COMM_WORLD,mystat,ierr)
            endif
            ! receive upper subdomain data and send lower subdomain data
            if (myid > 0) then
                call MPI_RECV(lower(:),nx,MPI_DOUBLE_PRECISION,myid-1,tag,MPI_COMM_WORLD,mystat,ierr)
                call MPI_SEND(Epot(:,1),nx,MPI_DOUBLE_PRECISION,myid-1,tag,MPI_COMM_WORLD,ierr)
            endif
        enddo

        ! get x-electric field from negative gradient of potential
        ! one-sided difference at left/right boundaries
        Ex(1,:) = -(Epot(2,:)-Epot(1,:))/dx
        Ex(nx,:) = -(Epot(nx,:)-Epot(nx-1,:))/dx
        ! central difference at interior points
        do i = 2, nx-1
            Ex(i,:) = -(Epot(i+1,:)-Epot(i-1,:))/(2*dx)
        enddo

        ! get y-electric field from negative gradient of potential
        do j = 1, nysub
            ! one-sided difference at lower boundary
            if (j == 1 .and. myid == 0) then
                Ey(:,1) = -(Epot(:,2)-Epot(:,1))/dy
            ! one-sided difference at upper boundary
            elseif (j == nysub .and. myid == nprc-1) then
                Ey(:,nysub) = -(Epot(:,nysub)-Epot(:,nysub-1))/dy
            ! central difference at interior points
            else
                if (j == 1) then
                    Ey(:,j) = -(Epot(:,j+1)-lower(:))/(2*dy)
                elseif (j == nysub) then
                    Ey(:,j) = -(upper(:)-Epot(:,j-1))/(2*dy)
                else
                    Ey(:,j) = -(Epot(:,j+1)-Epot(:,j-1))/(2*dy)
                endif
            endif
        enddo

        ! gather all E-field sub-domains across processes
        call MPI_ALLGATHER(Ex,nx*nysub,MPI_DOUBLE_PRECISION,allEx,nx*nysub,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        call MPI_ALLGATHER(Ey,nx*nysub,MPI_DOUBLE_PRECISION,allEy,nx*nysub,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        ! compute forces on particles from E-field on grid
        ! *Note: all processes see the E-field grid since particles migrate*
        do i = 1, myprt
            ndx = int(pxn(i)/dx)+1
            ndy = int(pyn(i)/dy)+1
            if (ndx < 1)  ndx = 1
            if (ndy < 1)  ndy = 1
            if (ndx > nx) ndx = nx
            if (ndy > ny) ndy = ny
            Fxn(i) = qn(i)*allEx(ndx,ndy)
            Fyn(i) = qn(i)*allEy(ndx,ndy)
        enddo

        ! advance particles using Verlet
        pxn(:) = pxn(:) + vxn(:)*dt + 0.5*Fxn(:)*dt*dt
        pyn(:) = pyn(:) + vyn(:)*dt + 0.5*Fyn(:)*dt*dt
        vxn(:) = vxn(:) + 0.5*Fxn(:)*dt
        vyn(:) = vyn(:) + 0.5*Fyn(:)*dt

        ! output first frame to file
        if (step == 1) then
            ! gather grid values to root array
            call MPI_GATHER(ygrid,nysub,MPI_DOUBLE_PRECISION,ygrid,nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(Epot,nx*nysub,MPI_DOUBLE_PRECISION,Epot,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(Ex,nx*nysub,MPI_DOUBLE_PRECISION,Ex,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(Ey,nx*nysub,MPI_DOUBLE_PRECISION,Ey,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            ! gather particles to root array
            call MPI_GATHER(pxn,myprt,MPI_DOUBLE_PRECISION,pxr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(pyn,myprt,MPI_DOUBLE_PRECISION,pyr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(vxn,myprt,MPI_DOUBLE_PRECISION,vxr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(vyn,myprt,MPI_DOUBLE_PRECISION,vyr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            ! only root process does this
            if (myid == 0) then
                ! write to file
                write(11,'(a8,f8.4)') 'time =',tstp
                write(12,'(a8,f8.4)') 'time =',tstp
                write(13,'(a8,f8.4)') 'time =',tstp
                write(14,'(5f12.6)') (tstp,pxr(i),pyr(i),vxr(i),vyr(i),i=1,nprt)
                do i = 1, nx
                    write(10,'(2f8.4)') xgrid(i),ygrid(i)
                    write(11,frmt) Epot(:,i)
                    write(12,frmt) Ex(:,i)
                    write(13,frmt) Ey(:,i)
                enddo
            endif
        endif

        ! output trajectory to file
        if (mod(step,nsav) == 0) then
            ! gather particles to root array
            call MPI_GATHER(pxn,myprt,MPI_DOUBLE_PRECISION,pxr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(pyn,myprt,MPI_DOUBLE_PRECISION,pyr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(vxn,myprt,MPI_DOUBLE_PRECISION,vxr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(vyn,myprt,MPI_DOUBLE_PRECISION,vyr,myprt,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            ! only root process does this
            if (myid == 0) then
                ! write to file
                write(14,'(5f12.6)') (tstp,pxr(i),pyr(i),vxr(i),vyr(i),i=1,nprt)
            endif
        endif
    enddo

    ! only root process does this
    if (myid == 0) then
        ! end timer
        call system_clock(t2,rate)
        ! print
        write(*,'(a25,i8,a5,i8)'),'grid =',nx,'x',ny
        write(*,'(a25,i8)'),'processors =',nprc
        write(*,'(a25,i8)'),'particles =',nprt
        write(*,'(a25,i8)'),'steps =',step-1
        write(*,'(a25,f8.3)'),'runtime =',real(t2-t1)/real(rate)
        ! deallocate arrays
        deallocate(pxr,pyr,vxr,vyr)
    endif

    ! deallocate arrays
    deallocate(xgrid,ygrid,chrg,Epot,Ex,Ey)
    deallocate(pxn,pyn,vxn,vyn,Fxn,Fyn,qn)
    deallocate(tmpchrg,allEx,allEy,upper,lower)

    ! end processes
    call MPI_FINALIZE(ierr)

end program particlemesh

! subroutine to generate positions
! generates random numbers over [0,1]
subroutine randnum(pxn,pyn,qn,myprt,myid)
    ! declare variables
    implicit none
    integer :: i,myprt,myid,t,clk
    real(8) :: a,xm,x0,seed
    real(8) :: pxn(myprt),pyn(myprt),qn(myprt)
    ! define constants
    a  = 16807.0
    xm = 2147483647.0
    x0 = 2147483711.0
    ! seed number generator
    call system_clock(t,clk)
    seed = clk*(myid+7)
    ! warmup number generator
    do i = 1, 200
        seed = mod(a*seed,xm)
    enddo
    ! start modulo generator
    do i = 1, myprt
        seed = mod(a*seed,xm)
        pxn(i) = seed/x0
    enddo
    do i = 1, myprt
        seed = mod(a*seed,xm)
        pyn(i) = seed/x0
    enddo
    do i = 1, myprt
        seed = mod(a*seed,xm)
        qn(i) = seed/x0
    enddo
end subroutine randnum