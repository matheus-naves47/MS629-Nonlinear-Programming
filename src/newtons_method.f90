module newtons_method

    implicit none

    abstract interface
        pure function fproc(x) result(result)
            double precision, dimension(:), intent(in) :: x ! Define como entrada o vetor x
            double precision                :: result ! Função retorna um valor
        end function fproc
    end interface

    abstract interface
        pure function gproc(x) result(result)
            double precision, dimension(:), intent(in) :: x ! Define como entrada o vetor x
            double precision, dimension(size(x)) :: result ! Função retorna um vetor
        end function gproc
    end interface

    abstract interface
        pure function hproc(x) result(result)
            double precision, dimension(:), intent(in) :: x ! Define como entrada o vetor x
            double precision, dimension(size(x), size(x)) :: result ! Função retorna uma matriz
        end function hproc
    end interface

contains

    subroutine metodo_newton(fun, grad, hess, x0, alpha, beta, gama, sigma, rho, eps, maxits)
        ! sub-rotina para otimizar uma função via método de Newton globalizado

        ! Declarações
        external :: dgesv ! Interface externa para sub-rotina DGESV do LAPACK
        procedure(fproc) :: fun ! Interface para função
        procedure(gproc) :: grad ! Interface para o gradiente da função
        procedure(hproc) :: hess ! Interface para a matriz hessiana da função
        double precision, dimension(:), intent(in) :: x0 ! Aproximação inicial x0
        double precision, dimension(size(x0), size(x0)) :: identidade

        double precision, intent(in) :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
        double precision, intent(in) :: beta ! Constante de proporcionalidade para a direção
        double precision, intent(in) :: gama ! Constante relativa ao ângulo entre a direção e o gradiente
        double precision, intent(in) :: sigma ! Fator de diminuição na busca unidimensional
        double precision, intent(in) :: rho ! Incremento inicial na globalização do método de Newton
        double precision, intent(in) :: eps ! Tolerância Epsilon

        integer, intent(in) :: maxits ! Número máximo de iterações
        integer :: i, j, k, n, info ! Contadores, tamanho do vetor x0 e status da resolução do sistema linear via sub-rotina DGESV
        integer :: pivot(size(x0)) ! Define vetor de permutações para fatoração LU da sub-rotina DGESV
        double precision, dimension(size(x0)) :: xk, dk ! Vetor iterador xk e vetor da direção do método
        double precision :: tk, mu

        info = 1
        k = 0
        xk = x0
        n = size(x0) ! Dimensão do vetor x0
        dk = 0

        ! Loop para definir a matriz identidade
        do i = 1, n
            do j = 1, n
                identidade(i, j) = 0
            end do
            identidade(i, i) = 1
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

        if (rho .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        if (eps .le. 0) then
            print *, "Parâmetros de entrada não aceitos, tente novamente!"
            stop
        end if

        do while (norm2(grad(xk)) .gt. eps) ! Critério de parada (precisão da norma do gradiente)
            mu = 0
            dk = -grad(xk)

            call dgesv(n, 1, (hess(xk) + mu*identidade), n, pivot, dk, n, info)

            if ((info .ne. 0) .or. (dot_product(grad(xk), dk) .gt. -gama*norm2(grad(xk))*norm2(dk))) then
                do while ((info .ne. 0) .or. (dot_product(grad(xk), dk) .gt. -gama*norm2(grad(xk))*norm2(dk)))
                    mu = max(2*mu, rho)
                    call dgesv(n, 1, (hess(xk) + mu*identidade), n, pivot, dk, n, info)
                end do
            end if

            if (norm2(dk) .lt. beta*norm2(grad(xk))) then
                dk = beta*(norm2(grad(xk))/norm2(dk))*dk
            end if

            tk = 1

            ! Verifica se o tamanho do passo respeita a regra de Armijo
            do while (fun(xk + tk*dk) .gt. (fun(xk) + alpha*tk*dot_product(grad(xk), dk)))
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
                    ! stop
                end if
            end if

        end do

        ! Impressão da saída
20      print 2, "******** Método de Newton ********"

        print *, "Parâmetros: "
        print 2, "alpha = ", alpha
        print 2, "beta = ", beta
        print 2, "gamma = ", gama
        print 2, "sigma = ", sigma
        print 2, "rho = ", rho
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
4       format(1X, "Encontrado após ", I0, " iterações.") ! Formatação do número de iterações.
5       format(1X, A, D24.17) ! Formatação do valor da função na solução

    end subroutine metodo_newton

end module newtons_method
