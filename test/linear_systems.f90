Program linear_system

    ! Declarações
    external :: dgesv ! Subrotina do LAPACK para resolução de sistemas lineares
    double precision :: A(3, 3), AF(3, 3), b(3), x(3), R(3), C(3) ! Define matriz a e o vetor b
    integer :: i, pivot(3), n, ok
    character(1) :: ch

    n = size(b)
    ch = "N"

    A(1, :) = (/1, 3, 5/)
    A(2, :) = (/2, 4, 6/)
    A(3, :) = (/1, 2, 7/)

    b(:) = (/4, 5, 7/)

    ! Teste com sistema sem solução
    ! A(1, :) = (/2, -1, 3/)
    ! A(2, :) = (/1, 2, -1/)
    ! A(3, :) = (/4, -2, 6/)

    ! b(:) = (/-4, 5, -8/)

    print *, b

    call dgesv(3, 1, A, 3, pivot, b, 3, ok)

    print *, ""

    print *, b
    print *, ok

    ! print the solution x
    do i = 1, 3
        write (*, 9) i, b(i)
    end do

9   format('x[', i1, '] = ', f6.3)
end program linear_system
