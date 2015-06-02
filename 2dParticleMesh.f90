! Created on Sat May 09 14:07:31 2015
! CSI 702, Assignment 6
! Nick Crump

! ** Serial version **

! 2D particle mesh simulation for electric charges
! Initializes particles to random positions and zero velocities
! Interpolates particle charge density to nearest grid points
! Solves Poisson equation using Jacobi/SOR iteration method
! Computes gradient of potential on grid to get E-field
! Interpolates forces on grid to nearest particles
! Advances particles using Verlet method
! Outputs values to text files

program particlemesh

    ! declare variables
    implicit none
    character            :: frmt*20,npts*7
    integer              :: i,j,step,nsav,nx,ny,iter,t1,t2,rate,nprt,nstp,ndx,ndy
    real(8)              :: pi,xmin,xmax,ymin,ymax,tol,tstp,dt
    real(8)              :: xi,yi,dx,dy,relax,mxerr,Eold,Enew
    real(8), allocatable :: xgrid(:),ygrid(:),chrg(:,:),Epot(:,:),Ex(:,:),Ey(:,:)
    real(8), allocatable :: pxn(:),pyn(:),vxn(:),vyn(:),Fxn(:),Fyn(:),qn(:)

    ! set input parameters
    ! -------------------------
    nprt = 200                ! number of particles
    nstp = 1000               ! number of time steps
    nsav = 10                 ! steps between saves
    tol  = 1e-5               ! accuracy tolerance
    xmin = 0    ; xmax = 20   ! x domain
    ymin = 0    ; ymax = 20   ! y domain
    nx   = 100  ;  ny  = 100  ! points in domain
    dt   = 0.005              ! time step
    ! -------------------------

    ! set parameters from input
    pi   = 4.0*atan(1.0)
    tstp = 0
    ! set relaxation parameter
    relax = 0.8
    ! set format write string for file output
    write(npts,'(i7)') nx
    frmt = '('//trim(npts)//'f12.6)'

    ! setup output files
    open(10,file='xygrid.txt')
    open(11,file='potential.txt')
    open(12,file='Exfield.txt')
    open(13,file='Eyfield.txt')
    open(14,file='trajectory.txt')

    ! start timer
    call system_clock(t1,rate)

    ! allocate arrays
    allocate(xgrid(nx),ygrid(ny),chrg(nx,ny),Epot(nx,ny),Ex(nx,ny),Ey(nx,ny))
    allocate(pxn(nprt),pyn(nprt),vxn(nprt),vyn(nprt),Fxn(nprt),Fyn(nprt),qn(nprt))

    ! get grid
    dx = (xmax-xmin)/(nx-1)
    dy = (ymax-ymin)/(ny-1)
    xi = xmin
    yi = ymin
    ! get x-grid
    do i = 1, nx
        xgrid(i) = xi
        xi = xi + dx
    enddo
    ! get y-grid
    do i = 1, ny
        ygrid(i) = yi
        yi = yi + dy
    enddo

    ! initialize particle positions and charges
    call randnum(pxn,pyn,qn,nprt)
    pxn(:) = (xmax-xmin)*pxn(:) + xmin
    pyn(:) = (ymax-ymin)*pyn(:) + ymin
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
        chrg(:,:)  = 0.0
        Epot(:,:)  = 0.0
        Ex(:,:)    = 0.0
        Ey(:,:)    = 0.0
        Fxn(:)     = 0.0
        Fyn(:)     = 0.0

        ! compute charge density on grid from charges on particles
        do i = 1, nprt
            ndx = int(pxn(i)/dx)+1
            ndy = int(pyn(i)/dy)+1
            if (ndx < 1)  ndx = 1
            if (ndy < 1)  ndy = 1
            if (ndx > nx) ndx = nx
            if (ndy > ny) ndy = ny
            chrg(ndx,ndy) = chrg(ndx,ndy) + qn(i)
        enddo

        ! iterate over potential until error is within tolerance
        do while (tol < mxerr)
            mxerr = 0
            iter  = iter + 1
            do j = 2, ny-1
                do i = 2, nx-1
                    ! get interior points of electric potential using SOR method
                    Eold = Epot(i,j)
                    Enew = 0.25*(Epot(i+1,j)+Epot(i-1,j)+Epot(i,j+1)+Epot(i,j-1)+chrg(i,j)*dx*dy)
                    Epot(i,j) = Enew + relax*(Enew-Eold)
                    ! accumulate instantaneous error
                    mxerr = mxerr + abs(Enew-Eold)
                enddo
            enddo
            ! get average error
            mxerr = mxerr/(nx*ny)
        enddo

        ! compute E-field on grid from gradient of potential
        do i = 2, nx-1
            ! central difference at xy interior points
            Ex(i,:) = -(Epot(i+1,:)-Epot(i-1,:))/(2*dx)
            Ey(:,i) = -(Epot(:,i+1)-Epot(:,i-1))/(2*dy)
        enddo

        ! compute forces on particles from E-field on grid
        do i = 1, nprt
            ndx = int(pxn(i)/dx)+1
            ndy = int(pyn(i)/dy)+1
            if (ndx < 1)  ndx = 1
            if (ndy < 1)  ndy = 1
            if (ndx > nx) ndx = nx
            if (ndy > ny) ndy = ny
            Fxn(i) = qn(i)*Ex(ndx,ndy)
            Fyn(i) = qn(i)*Ey(ndx,ndy)
        enddo

        ! advance particles using Verlet
        pxn(:) = pxn(:) + vxn(:)*dt + 0.5*Fxn(:)*dt*dt
        pyn(:) = pyn(:) + vyn(:)*dt + 0.5*Fyn(:)*dt*dt
        vxn(:) = vxn(:) + 0.5*Fxn(:)*dt
        vyn(:) = vyn(:) + 0.5*Fyn(:)*dt

        ! output first frame to file
        if (step == 1) then
            write(11,'(a8,f8.4)') 'time =',tstp
            write(12,'(a8,f8.4)') 'time =',tstp
            write(13,'(a8,f8.4)') 'time =',tstp
            write(14,'(5f12.6)') (tstp,pxn(i),pyn(i),vxn(i),vyn(i),i=1,nprt)
            do i = 1, nx
                write(10,'(2f8.4)') xgrid(i),ygrid(i)
                write(11,frmt) Epot(:,i)
                write(12,frmt) Ex(:,i)
                write(13,frmt) Ey(:,i)
            enddo
        endif

        ! output trajectory to file
        if (mod(step,nsav) == 0) then
            write(14,'(5f12.6)') (tstp,pxn(i),pyn(i),vxn(i),vyn(i),i=1,nprt)
        endif
    enddo

    ! end timer
    call system_clock(t2,rate)

    ! print
    write(*,'(a25,i8,a5,i8)'),'grid =',nx,'x',ny
    write(*,'(a25,i8)'),'particles =',nprt
    write(*,'(a25,i8)'),'steps =',step-1
    write(*,'(a25,f8.3)'),'runtime =',real(t2-t1)/real(rate)

    ! deallocate arrays
    deallocate(xgrid,ygrid,chrg,Epot,Ex,Ey)
    deallocate(pxn,pyn,vxn,vyn,Fxn,Fyn,qn)

end program particlemesh

! subroutine to generate positions
! generates random numbers over [0,1]
subroutine randnum(pxn,pyn,qn,nprt)
    ! declare variables
    implicit none
    integer :: i,nprt,clk
    real(8) :: a,xm,x0,seed
    real(8) :: pxn(nprt),pyn(nprt),qn(nprt)
    ! define constants
    a  = 16807.0
    xm = 2147483647.0
    x0 = 2147483711.0
    ! seed number generator
    call system_clock(count=clk)
    seed = clk+clk*1.2
    ! start modulo generator
    do i = 1, nprt
        seed = mod(a*seed,xm)
        pxn(i) = seed/x0
    enddo
    do i = 1, nprt
        seed = mod(a*seed,xm)
        pyn(i) = seed/x0
    enddo
    do i = 1, nprt
        seed = mod(a*seed,xm)
        qn(i) = seed/x0
    enddo
end subroutine randnum