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

end module functions
