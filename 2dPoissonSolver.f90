! Created on Thu Apr 30 15:06:23 2015
! CSI 702, Assignment 5
! Nick Crump

! ** Serial version **

! Solves Poisson equation for electric potential on square 2D grid
! Uses Jacobi/SOR iterative method to solve PDE
! Uses Dirichlet boundary conditions
! Calculates gradient of the solution to get electric field
! Outputs values to text files

program poissonsolve

    ! declare variables
    implicit none
    character            :: frmt*20,npts*7
    integer              :: i,j,nx,ny,iter,t1,t2,rate,p
    real(8)              :: pi,xmin,xmax,ymin,ymax,bx0,bx1,tol
    real(8)              :: xi,yi,dx,dy,k,relax,errs,Eold,Enew
    real(8), allocatable :: x(:),y(:),Epot(:,:),chrg(:,:),gradx(:,:),grady(:,:)

    ! set input parameters
    ! -------------------------
    tol  =  1e-5              ! accuracy tolerance
    xmin = -10  ; xmax = 10   ! x domain
    ymin = -10  ; ymax = 10   ! y domain
    nx   =  100 ;  ny  = 100  ! points in domains
    bx0  = -1.0 ; bx1  = 1.0  ! left/right boundary values
    ! -------------------------

    ! setup output files
    open(10,file='xygrid.txt')
    open(11,file='potential.txt')
    open(12,file='gradx.txt')
    open(13,file='grady.txt')

    ! start timer
    call system_clock(t1,rate)

    ! allocate arrays
    allocate(x(nx),y(ny),Epot(nx,ny),chrg(nx,ny),gradx(nx,ny),grady(nx,ny))

    ! initialize arrays to zero
    chrg(:,:)  = 0.0
    Epot(:,:)  = 0.0
    gradx(:,:) = 0.0
    grady(:,:) = 0.0

    ! set some parameters
    pi   = 4.0*atan(1.0)
    iter = 0
    errs = 1
    ! set format write string for file output
    write(npts,'(i7)') nx
    frmt = '('//trim(npts)//'f12.6)'
    ! set relaxation parameter
    relax = 0.8

    ! get grid
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)
    xi = xmin
    yi = ymin
    ! get x-grid
    do i = 1, nx
        x(i) = xi
        xi = xi + dx
    enddo
    ! get y-grid
    do i = 1, ny
        y(i) = yi
        yi = yi + dy
    enddo

    ! place single point charge at center of grid (RHS of Poisson equation)
    chrg(nx/2,ny/2) = 1.0/dx**2

    ! set left/right boundary values
    Epot(1,:)  = bx0
    Epot(nx,:) = bx1

    ! iterate over solution until error is within tolerance
    do while (tol < errs)
!        print*,'iter =',iter
!        write(*,frmt) (Epot(:,p),p=1,nx)
!        print*,''
        errs = 0
        iter = iter + 1
        do j = 1, ny
            do i = 2, nx-1
                ! get top/bottom boundary values
                if (j == 1 .or. j == ny) then
                    k = (x(i)-xmin)/(xmax-xmin)
                    Epot(i,j) = (1-k)*bx0 + k*bx1
                else
                    ! get interior points of electric potential using SOR method
                    Eold = Epot(i,j)
                    Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+Epot(i,j-1)+chrg(i,j)*dx**2)
                    Epot(i,j) = Enew + relax*(Enew-Eold)
                    ! accumulate instantaneous error
                    errs = errs + abs(Enew-Eold)
                endif
            enddo
        enddo
        ! get average error
        errs = errs/(nx*ny)
    enddo

    ! get electric field from negative gradient of potential
    ! forward difference at xy lower boundaries
    gradx(1,:) = -(Epot(2,:)-Epot(1,:))/dx
    grady(:,1) = -(Epot(:,2)-Epot(:,1))/dy
    ! backward difference at xy upper boundaries
    gradx(nx,:) = -(Epot(nx,:)-Epot(nx-1,:))/dx
    grady(:,ny) = -(Epot(:,ny)-Epot(:,ny-1))/dy
    ! central difference at xy interior points
    do i = 2, ny-1
        gradx(i,:) = -(Epot(i+1,:)-Epot(i-1,:))/(2*dx)
        grady(:,i) = -(Epot(:,i+1)-Epot(:,i-1))/(2*dy)
    enddo

    ! end timer
    call system_clock(t2,rate)

    ! print
    write(*,'(a25,i8,a5,i8)'),'grid size =',nx,'x',ny
    write(*,'(a25,i8)'),'iterations to converge =',iter
    write(*,'(a25,f8.3)'),'seconds to converge =',real(t2-t1)/real(rate)

    ! output to file
    do i = 1, ny
        write(10,*) x(i),y(i)
        write(11,frmt) Epot(:,i)
        write(12,frmt) gradx(:,i)
        write(13,frmt) grady(:,i)
    enddo

    ! deallocate arrays
    deallocate(x,y,Epot,chrg,gradx,grady)

end program poissonsolve