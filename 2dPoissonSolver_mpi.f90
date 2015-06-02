! Created on Thu Apr 30 15:06:23 2015
! CSI 702, Assignment 5
! Nick Crump

! ** MPI version **

! Solves Poisson equation for electric potential on square 2D grid
! Subdivides domain along y-axis into blocks for each process
! Uses Jacobi/SOR iterative method to solve PDE
! Uses Dirichlet boundary conditions
! Calculates gradient of the solution to get electric field
! Outputs values to text files

program poissonsolve

    ! import mpi
    use mpi
    ! declare variables
    implicit none
    character            :: frmt*20,npts*7
    integer              :: ierr,nprc,myid,tag,mystat(MPI_STATUS_SIZE)
    integer              :: i,j,p,nx,ny,nysub,iter,mxiter,t1,t2,rate
    real(8)              :: pi,xmin,xmax,ymin,ymax,bx0,bx1,tol
    real(8)              :: dx,dy,k,relax,myerr,mxerr,Eold,Enew,Lx,Ly
    real(8), allocatable :: x(:),y(:),chrg(:,:),Epot(:,:),Ex(:,:),Ey(:,:)
    real(8), allocatable :: upper(:),lower(:)

    ! set input parameters
    ! ------------------------------------------------------
    tol    =  1e-5                ! accuracy tolerance
    mxiter = 10000                ! max iterations
    xmin   = -10   ; xmax = 10    ! x domain
    ymin   = -10   ; ymax = 10    ! y domain
    nx     =  2000 ;  ny  = 2000  ! points in domains
    bx0    = -1.0  ; bx1  = 1.0   ! left/right boundary values
    ! ------------------------------------------------------

    ! set parameters from input
    pi    = 4.0*atan(1.0)
    iter  = 0
    mxerr = 1
    Lx = xmax-xmin
    Ly = ymax-ymin
    ! set relaxation parameter
    relax = 0.8

    ! get mpi processes and ids
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprc,ierr)

    ! get domain partition (y-axis subdivided)
    nysub = ny/nprc
    dx    = Lx/(nx-1)
    dy    = Ly/(ny-1)

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
        open(12,file='gradx.txt')
        open(13,file='grady.txt')
        ! set format write string for output files
        write(npts,'(i7)') nx
        frmt = '('//trim(npts)//'f12.6)'
        ! start timer
        call system_clock(t1,rate)
        ! allocate domain arrays to root process
        allocate(x(nx),y(ny),chrg(nx,ny),Epot(nx,ny),Ex(nx,ny),Ey(nx,ny),upper(nx),lower(nx))
    ! all other processes do this
    else
       ! allocate subdomain arrays to all others
       allocate(x(nx),y(nysub),chrg(nx,nysub),Epot(nx,nysub),Ex(nx,nysub),Ey(nx,nysub),upper(nx),lower(nx))
    endif

    ! initialize arrays to zero
    chrg(:,:)  = 0.0
    Epot(:,:)  = 0.0
    Ex(:,:)    = 0.0
    Ey(:,:)    = 0.0
    upper(:)   = 0.0
    lower(:)   = 0.0

    ! get x-grid
    do i = 1, nx
        x(i) = xmin + dx*(i-1)
    enddo
    ! get y-grid (subdivided)
    do j = 1, nysub
        y(j) = ymin + dy*(j-1) + dy*myid*nysub
    enddo

    ! place single point charge at center of grid
    ! only process owning center point does this
    if (myid == nprc/2) then
        chrg(nx/2,1) = 1.0/dx**2
    endif

    ! set left/right boundary values
    Epot(1,:)  = bx0
    Epot(nx,:) = bx1

    ! iterate over solution until error is within tolerance
    do while (tol < mxerr)
        ! stop if max iterations reached
        if (iter == mxiter) then
            if (myid == 0) then
                print*, 'max iterations reached =',iter
                call MPI_FINALIZE(ierr)
                stop
            endif
        endif
        myerr = 0
        iter  = iter + 1
        do j = 1, nysub
            do i = 2, nx-1
                ! get bottom boundary value
                if (j == 1 .and. myid == 0) then
                    k = (x(i)-xmin)/Lx
                    Epot(i,j) = (1-k)*bx0 + k*bx1
                ! get top boundary value
                elseif (j == nysub .and. myid == nprc-1) then
                    k = (x(i)-xmin)/Lx
                    Epot(i,j) = (1-k)*bx0 + k*bx1
                else
                    Eold = Epot(i,j)
                    ! get interior points of electric potential using SOR method
                    if (j == 1) then
                        Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+lower(i)+chrg(i,j)*dx**2)
                    elseif (j == nysub) then
                        Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+upper(i)+Epot(i,j-1)+chrg(i,j)*dx**2)
                    else
                        Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+Epot(i,j-1)+chrg(i,j)*dx**2)
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
        call MPI_Allreduce(myerr,mxerr,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
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

    ! gather all values to root array
    call MPI_GATHER(y,nysub,MPI_DOUBLE_PRECISION,y,nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(Epot,nx*nysub,MPI_DOUBLE_PRECISION,Epot,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(Ex,nx*nysub,MPI_DOUBLE_PRECISION,Ex,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(Ey,nx*nysub,MPI_DOUBLE_PRECISION,Ey,nx*nysub,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! only root process does this
    if (myid == 0) then
        ! end timer
        call system_clock(t2,rate)
        ! print
        write(*,'(a25,i8,a5,i8)'),'grid size =',nx,'x',ny
        write(*,'(a25,i8)'),'number of processors =',nprc
        write(*,'(a25,i8)'),'iterations to converge =',iter
        write(*,'(a25,f8.3)'),'seconds to converge =',real(t2-t1)/real(rate)
        ! output to file
        do i = 1, ny
            write(10,*) x(i),y(i)
            write(11,frmt) Epot(:,i)
            write(12,frmt) Ex(:,i)
            write(13,frmt) Ey(:,i)
        enddo
    endif

    ! deallocate arrays
    deallocate(x,y,chrg,Epot,Ex,Ey,upper,lower)
    ! end processes
    call MPI_FINALIZE(ierr)

end program poissonsolve