program hessian_approx
    !
    use gradient_method
    use functions
    implicit none

    real :: hessian(4, 4)
    real, dimension(4) :: z
    ! integer, parameter::n = size(hessian)
    ! integer, parameter::n = shape(hessian)
    integer :: i, j, n

    n = size(z)

    ! z = (/1., 2., 3./)

    print *, "Teste Matriz Hessiana (Aproximação)!"

    do i = 1, n
        do j = 1, n
            hessian(i, j) = 0.
        end do
        hessian(i, i) = 1.
    end do

    do i = 1, n
        print *, hessian(i, :)
    end do

end program hessian_approx
