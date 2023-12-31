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
        double precision, dimension(size(x0)) :: xk, dk ! Vetor iterador xk e vetor da direção do método
        double precision :: tk

        k = 0 ! Inicializa contador de iterações
        xk = x0 ! Inicializa xk como x0
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
            dk = -grad(xk) ! Toma direção de descida como gradiente conjugado
            tk = 1. ! Iterador para regra de Armijo

            ! Verifica se o tamanho do passo respeita a regra de Armijo
            do while (fun(xk + tk*dk) .gt. (fun(xk) + alpha*tk*dot_product(dk, -dk)))
                tk = sigma*tk ! Tamanho do passo
            end do
            xk = xk + tk*dk ! Atualiza o vetor xk
            k = k + 1 ! Atualiza o contador de iterações

            ! Se maxits for definido como 0
            ! o método é executado apenas
            ! com o critério de parada do epsilon.
            if (maxits .ne. 0) then
                if (k == maxits) then
                    goto 20
                end if
            end if

        end do

        ! Impressão da saída
20      print 2, "******** Método do Gradiente ********"

        print *, "Parâmetros: "
        print 2, "alpha = ", alpha
        print 2, "sigma = ", sigma
        print 2, "eps = ", eps
        print *, ""
        print *, "Solução encontrada em:"
        ! Imprime cada valor de xi em cada linha
        do i = 1, n
            write (*, 3) i, xk(i)
        end do
        print *, ""
        write (*, 5) "Função avaliada na solução: ", fun(xk)
        write (*, 4) k
2       format(1X, A, D0.2) ! Formatação dos parâmetros
3       format(1X, 'x[', I0, '] = ', D24.17) ! Formatação de cada xi
4       format(1X, "Encontrada após ", I0, " iterações.") ! Formatação do número de iterações.
5       format(1X, A, D24.17) ! Formatação do valor da função na solução

    end subroutine metodo_gradiente
end module gradient_method
