program ompArray
    implicit none
    integer :: i,n = 1000
    real    :: s(1000)
    !$omp parallel do private(i) shared(n,s)
    do i = 1,n
        s(i) = i
    enddo
    !$omp end parallel do
    print*, sum(s)
end program ompArray