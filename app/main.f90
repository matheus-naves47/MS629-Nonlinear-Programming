program OTISER
    ! Programa OTImização SEm Restrições
    use gradient_method ! Módulo do método do Gradiente
    use newtons_method ! Módulo do método de Newton Globalizado
    use cp1_method ! Módulo do método de Correção de Posto 1
    use dfp_method ! Módulo do método de Davidon, Fletcher e Powell
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    external :: DGESV ! Subrotina do LAPACK para resolver sistema linear
    double precision, dimension(6):: x, y, z, w
    integer, parameter::n = size(x)
    double precision :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
    double precision :: beta ! Constante de proporcionalidade para a direção
    double precision :: gama ! Constante relativa ao ângulo entre a direção e o gradiente
    double precision :: sigma ! Fator de diminuição na busca unidimensional
    double precision :: rho ! Incremento inicial na globalização do método de Newton
    double precision :: eps ! precisão para a norma do gradiente (critério de parada).
    integer :: i, maxits ! Contador e número máximo de iterações

    alpha = 1e-4
    beta = 1e-3
    gama = 1e-6
    sigma = 0.5
    rho = 1e-3
    eps = 1e-2
    maxits = 479760

    ! Definição do vetor inicial x0
    do i = 1, n
        x(i) = 1
        y(i) = 0
        z(i) = i
        w(i) = i*i
    end do

    ! Chama os métodos na função quadrática. Pode-se testar com outras funções e outros parâmetros. 
    ! Método do Gradiente
    call metodo_gradiente(quadratica, quadratica_grad, x, alpha, sigma, eps, maxits)
    print *, ""

    ! Método de Newton
    call metodo_newton(quadratica, quadratica_grad, quadratica_hess, x, alpha, beta, gama, sigma, rho, eps, maxits)
    print *, ""

    ! Método CP1
    call metodo_cp1(quadratica, quadratica_grad, x, alpha, beta, gama, sigma, eps, maxits)
    print *, ""

    ! Método DFP
    call metodo_dfp(quadratica, quadratica_grad, x, alpha, beta, gama, sigma, eps, maxits)
    print *, ""

end program OTISER
