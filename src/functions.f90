module functions

contains

    ! Exemplo 1 - Quadrática
    pure function quadratica(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result = result + (i*(x(i)**2))
        end do
    end function quadratica

    pure function quadratica_grad(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x)) :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result(i) = (2*i*x(i))
        end do
    end function quadratica_grad

    pure function quadratica_hess(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x), size(x)) :: result
        integer :: i, j, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0

        do i = 1, n ! loop para definir a matriz hessiana
            do j = 1, n
                result(i, j) = 0
            end do
            result(i, i) = 2*i
        end do
    end function quadratica_hess

    ! Exemplo 2 - Rosenbrock
    pure function rosenbrock(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, (n/2) ! loop para definir a função
            result = result + (100*(x(2*i) - (x(2*i - 1))**2)**2 + (x(2*i - 1) - 1)**2)
        end do
    end function rosenbrock

    pure function rosenbrock_grad(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x)) :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            if (mod(i, 2) .ne. 0) then
                result(i) = 2*(x(i) - 1) - 400*x(i)*(x(i + 1) - x(i)**2)
            else if (mod(i, 2) .eq. 0) then
                result(i) = 200*(x(i) - x(i - 1)**2)
            end if
        end do
    end function rosenbrock_grad

    pure function rosenbrock_hess(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x), size(x)) :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            if (mod(i, 2) .ne. 0) then
                result(i,i) = -20
            else if (mod(i, 2) .eq. 0) then
                result(i,i) = 2
            end if
        end do
    end function rosenbrock_hess

    ! Exemplo 3 - Styblinsky-Tang
    pure function styblinsky_tang(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result = result + ((x(i)*x(i)*x(i)*x(i)) - 16*(x(i)*x(i)) + 5*x(i))
        end do
    end function styblinsky_tang

    pure function styblinsky_tang_grad(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x)) :: result
        integer :: i, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result(i) = (4*(x(i)*x(i)*x(i)) - 32*(x(i)) + 5)
        end do
    end function styblinsky_tang_grad

    pure function styblinsky_tang_hess(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x), size(x)) :: result
        integer :: i, j, n ! índices do somatório
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0

        do i = 1, n ! loop para definir a matriz hessiana
            do j = 1, n
                result(i, j) = 0
            end do
            result(i, i) = (12*(x(i)*x(i)) - 32)
        end do
    end function styblinsky_tang_hess

    ! Exemplo 4 - Rastrigin
    pure function rastrigin(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision :: result
        integer :: i, n ! índices do somatório
        double precision :: pi ! Define Pi em precisão dupla
        pi = 4*atan(1.0d0) ! Fórmula da tangente para Pi em precisão dupla
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result = result + (x(i)*(x(i)) - 10*cos(2*pi*x(i)))
        end do
    end function rastrigin

    pure function rastrigin_grad(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x)) :: result
        integer :: i, n ! índices do somatório
        double precision :: pi ! Define Pi em precisão dupla
        pi = 4*atan(1.0d0) ! Fórmula da tangente para Pi em precisão dupla
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0
        do i = 1, n ! loop para definir a função
            result(i) = 2*x(i) + 20*pi*sin(2*pi*x(i))
        end do
    end function rastrigin_grad

    pure function rastrigin_hess(x) result(result)

        ! Declarações
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x), size(x)) :: result
        integer :: i, j, n ! índices do somatório
        double precision :: pi ! Define Pi em precisão dupla
        pi = 4*atan(1.0d0)
        n = size(x) ! n é o tamanho do vetor de entrada
        result = 0 ! Função inicializa no 0

        do i = 1, n ! loop para definir a matriz hessiana
            do j = 1, n
                result(i, j) = 0
            end do
            result(i, i) = 2 + 40*pi*pi*cos(2*pi*x(i))
        end do

    end function rastrigin_hess

end module functions
