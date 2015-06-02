program ompSum
    implicit none
    integer :: i,n,s
    s = 0
    n = 100000000
    !$omp parallel do private(i) shared(n) reduction(+:s)
    do i = 1,n
        s = s + i
    enddo
    !$omp end parallel do
    print*, s
end program ompSum