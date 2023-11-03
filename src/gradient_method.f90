module gradient_method
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

contains

    subroutine metodo_gradiente(fun, grad, x0, alpha, sigma, eps, maxits)
        ! Subrotina para otimizar uma função via método do gradiente

        procedure(fproc) :: fun
        procedure(gproc) :: grad
        double precision, dimension(:), intent(in) :: x0 ! Aproximação inicial x0
        double precision, intent(in) :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
        double precision, intent(in) :: sigma ! Fator de diminuição na busca unidimensional
        double precision, intent(in) :: eps ! Tolerância Epsilon
        integer, intent(in) :: maxits ! Número máximo de iterações
        integer :: i, k, n ! Contadores e tamanho do vetor x0
        double precision, dimension(size(x0)) :: xk, direcao ! Vetor iterador xk e vetor da direção do método
        double precision :: tk

        k = 0
        xk = x0
        n = size(x0) ! Define n como a dimensão do vetor x0

        ! Verificação de parâmetros

        if (alpha .ge. 1 .or. alpha .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (sigma .ge. 1 .or. sigma .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (eps .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        do while (norm2(grad(xk)) .gt. eps) ! Critério de parada (precisão da norma do gradiente)
            direcao = -grad(xk) ! Toma direção de descida como gradiente conjugado
            tk = 1. ! Iterador para regra de Armijo

            ! Verifica se o tamanho do passo respeita a regra de Armijo
            do while (fun(xk + tk*direcao) .gt. (fun(xk) + alpha*tk*dot_product(direcao, -direcao)))
                tk = sigma*tk ! Tamanho do passo
            end do
            xk = xk + tk*direcao ! Atualiza o vetor xk
            k = k + 1 ! Atualiza o contador de iterações

            ! Se maxits for definido como 0
            ! o método é executado apenas
            ! com o critério de parada do epsilon.
            if (maxits .ne. 0) then
                if (k == maxits) then
                    goto 10
                end if
            end if

        end do

! 10      print 2, "Mínimo local em: ", ("X", i, "=", xk(i), i=1, size(xk))
! 2       format(1X, A, 1X, (5(1X, A, I0, A, 1X, F7.2)))
10      print 2, "******** Método do Gradiente ********"

        print *, "Parâmetros: "
        print 2, "alpha = ", alpha
        print 2, "sigma = ", sigma
        print 2, "epsilon = ", eps
2       format(1X, A, ES10.2E2)

        print *, ""

        print *, "Mínimo local em:"
        do i = 1, n
            print 3, "x", i, " = ", xk(i)
3           format(1X, A, I0, A, F22.16)
! 2           format(1X, A, I0, A, ES15.3E3)
        end do
        print 4, "Encontrado após ", k, " iterações."
4       format(1X, A, I12, A)

        print *, ""

    end subroutine metodo_gradiente
end module gradient_method
