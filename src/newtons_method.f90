module newtons_method
    abstract interface
        pure function fproc(x) result(result)
            double precision, dimension(:), intent(in) :: x
            double precision                :: result
        end function fproc
    end interface

    abstract interface
        pure function gproc(x) result(result)
            double precision, dimension(:), intent(in) :: x
            double precision, dimension(size(x)) :: result
        end function gproc
    end interface

    abstract interface
        pure function hproc(x) result(result)
            double precision, dimension(:), intent(in) :: x
            double precision, dimension(size(x), size(x)) :: result
        end function hproc
    end interface

contains

    subroutine metodo_newton(fun, grad, hess, x0, alpha, beta, gama, sigma, rho, eps, maxits)
        ! subrotina para otimizar uma função via método de Newton globalizado

        ! Declarações
        procedure(fproc) :: fun
        procedure(gproc) :: grad
        procedure(hproc) :: hess
        double precision, dimension(:), intent(in) :: x0 ! Aproximação inicial x0
        double precision, intent(in) :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
        double precision :: beta ! Constante de proporcionalidade para a direção
        double precision :: gama ! Constante relativa ao ângulo entre a direção e o gradiente
        double precision, intent(in) :: sigma ! Fator de diminuição na busca unidimensional
        double precision :: rho ! Incremento inicial na globalização do método de Newton
        double precision, intent(in) :: eps ! Tolerância Epsilon
        integer, intent(in) :: maxits ! Número máximo de iterações
        integer :: i, k, n ! Contadores e tamanho do vetor x0
        double precision, dimension(size(x0)) :: xk, direcao ! Vetor iterador xk e vetor da direção do método
        double precision :: tk, mu

        k = 0
        xk = x0
        n = size(x0) ! Dimensão do vetor x0

        ! Verificaçao de parâmetros

        if (alpha .ge. 1 .or. alpha .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (beta .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (gama .ge. 1 .or. gama .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (sigma .ge. 1 .or. sigma .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (rho .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (eps .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        do while (norm2(grad(fun, xk)) .gt. eps) ! Critério de parada (precisão da norma do gradiente)
            mu = 0
        end do

    end subroutine metodo_newton

end module newtons_method
