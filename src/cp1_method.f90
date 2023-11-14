module cp1_method
    implicit none

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

    subroutine metodo_cp1(fun, grad, x0, alpha, beta, gama, sigma, eps, maxits)
        ! Subrotina para otimizar uma função via método quasi-newton Correção de Posto 1

        ! Declarações
        procedure(fproc) :: fun ! Interface para função
        procedure(gproc) :: grad ! Interface para o gradiente da função
        double precision, dimension(:), intent(in) :: x0 ! Aproximação inicial x0
        double precision, dimension(size(x0), size(x0)) :: hk ! Aproximação inicial da inversa da hessiana (identidade)

        double precision, intent(in) :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
        double precision, intent(in) :: beta ! Constante de proporcionalidade para a direção
        double precision, intent(in) :: gama ! Constante relativa ao ângulo entre a direção e o gradiente
        double precision, intent(in) :: sigma ! Fator de diminuição na busca unidimensional
        double precision, intent(in) :: eps ! Tolerância Epsilon

        integer, intent(in) :: maxits ! Número máximo de iterações
        integer :: i, j, k, n ! Contadores e dimensão do vetor x0
        double precision, dimension(size(x0)) :: xk, dk ! Vetor iterador xk e vetor da direção do método
        double precision, dimension(size(x0)) :: sk, yk, zk, wk ! Vetores auxiliares
        double precision :: tk ! Iterador para regra de armijo

        ! Definições
        k = 0
        xk = x0
        n = size(x0)

        ! Loop para definir a matriz identidade
        do i = 1, n
            do j = 1, n
                hk(i, j) = 0
            end do
            hk(i, i) = 1
        end do

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

        if (eps .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        do while (norm2(grad(xk)) .gt. eps) ! Critério de parada (precisão da norma do gradiente)

            dk = matmul(-hk, grad(xk))

            if (dot_product(grad(xk), dk) .gt. -gama*norm2(grad(xk))*norm2(dk)) then
                dk = -grad(xk)

                do i = 1, n
                    do j = 1, n
                        hk(i, j) = 0
                    end do
                    hk(i, i) = 1
                end do

            end if

            if (norm2(dk) .lt. beta*norm2(grad(xk))) then
                dk = beta*(norm2(grad(xk))/norm2(dk))*dk
            end if

            tk = 1

            do while (fun(xk + tk*dk) .gt. (fun(xk) + alpha*tk*dot_product(grad(xk), dk)))
                tk = sigma*tk ! Tamanho do passo
            end do

            ! Define vetores auxiliares e atualiza o vetor xk
            sk = (xk + tk*dk) - xk
            yk = grad(xk + tk*dk) - grad(xk)
            zk = matmul(hk, yk)
            wk = sk - zk
            xk = xk + tk*dk

            ! Atualiza a matriz hk

            if (dot_product(wk, yk) .gt. 0) then
                hk = hk + (dot_product(wk, wk)/dot_product(wk, yk))
                ! print *, hk
            else
                hk = hk
                ! print *, hk
            end if

            k = k + 1 ! Atualiza o contador de iterações

            ! Se maxits for definido como 0
            ! o método é executado apenas
            ! com o critério de parada do epsilon.
            if (maxits .ne. 0) then
                if (k == maxits) then
                    goto 20
                    ! stop
                end if
            end if

        end do

20      print 2, "******** Método de Correção de Posto 1 ********"

        print *, "Parâmetros: "
        print 2, "alpha = ", alpha
        print 2, "beta = ", beta
        print 2, "gamma = ", gama
        print 2, "sigma = ", sigma
        print 2, "epsilon = ", eps

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
    end subroutine metodo_cp1

end module cp1_method
