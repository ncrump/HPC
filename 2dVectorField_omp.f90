! Created on Sat Feb 21 20:09:00 2015
! CSI 702, Assignment 3
! Nick Crump

! ** OpenMP version**

! Calculates gradient of 2D scalar field using finite differences
! Outputs values to text files
program vectorfield
    ! declare variables
    implicit none
    integer,parameter   :: n = 99
    integer             :: i
    real                :: gmin,gmax,h,pxy,magmax
    real,dimension(n)   :: xgrid,ygrid
    real,dimension(n,n) :: ff,gx,gy,mag
!    ! setup output files
!    open(10,file='xygrid.txt')
!    open(11,file='2dscalarfield.txt')
!    open(12,file='2dvectorfieldX.txt')
!    open(13,file='2dvectorfieldY.txt')
    ! set limits
    gmin = -1.0 ; gmax = 1.0
    ! set arrays to zero
    ff(:,:) = 0.0
    gx(:,:) = 0.0
    gy(:,:) = 0.0
    ! get grid
    h = (gmax-gmin)/n
    pxy = gmin
    do i = 1, n
        xgrid(i) = pxy
        ygrid(i) = pxy
        pxy = pxy + h
    enddo

    ! get functions
    call func2d(1.0,1.0,n,xgrid,ygrid,ff)
    call grad2d(n,h,ff,gx,gy)
    call maxmag(n,gx,gy,mag,magmax)
    ! print values
    print*,'max magnitude of gradient =',magmax
!    ! output to file
!    do i = 1,n
!    write(10,*) xgrid(i),ygrid(i)
!    write(11,*) ff(:,i)
!    write(12,*) gx(:,i)
!    write(13,*) gy(:,i)
!    enddo
end program vectorfield

! subroutine to calculate 2d scalar field
subroutine func2d(a,b,n,xgrid,ygrid,ff)
    ! declare variables
    implicit none
    integer :: i,j,n
    real    :: a,b
    real    :: xgrid(n),ygrid(n),ff(n,n)
    ! loop to calculate function
    !$omp parallel do private(i,j) shared(a,b,n,xgrid,ygrid,ff)
    do i = 1, n
        do j = 1, n
            ff(i,j) = a*xgrid(i) + (b*xgrid(i))/(xgrid(i)**2 + ygrid(j)**2)
        enddo
    enddo
    !$omp end parallel do
end subroutine func2d

! subroutine to calculate 2d gradient of scalar field
subroutine grad2d(n,h,ff,gx,gy)
    ! declare variables
    integer :: i,n
    real    :: h
    real    :: ff(n,n),gx(n,n),gy(n,n)
    ! forward difference at xy lower boundaries
    gx(:,1) = (ff(:,2)-ff(:,1))/h
    gy(1,:) = (ff(2,:)-ff(1,:))/h
    ! backward difference at xy upper boundaries
    gx(:,n) = (ff(:,n)-ff(:,n-1))/h
    gy(n,:) = (ff(n,:)-ff(n-1,:))/h
    ! central difference at xy interior points
    !$omp parallel do private(i) shared(n,h,gx,gy,ff)
    do i = 2, n-1
        gx(:,i) = (ff(:,i+1)-ff(:,i-1))/(2*h)
        gy(i,:) = (ff(i+1,:)-ff(i-1,:))/(2*h)
    enddo
    !$omp end parallel do
end subroutine grad2d

! subroutine to calculate maximum magnitude of gradient
subroutine maxmag(n,gx,gy,mag,magmax)
    ! declare variables
    implicit none
    integer :: i,j,n
    real    :: magmax
    real    :: gx(n,n),gy(n,n),mag(n,n)
    ! loop to calculate magnitude
    !$omp parallel do private(i) shared(n,gx,gy,mag)
    do i = 1, n
        mag(i,:) = (gx(i,:)**2 + gy(i,:)**2)**0.5
    enddo
    !$omp end parallel do
    ! loop to calculate max
    magmax = 0
    !$omp parallel do private(i,j) shared(n,mag,magmax)
    do i = 1, n
        do j = 1, n
            if (mag(i,j) > magmax) then
                magmax = mag(i,j)
            endif
        enddo
    enddo
    !$omp end parallel do
end subroutine maxmag