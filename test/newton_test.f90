program newton
    use newtons_method ! Módulo do método de Newton Globalizado
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    ! double precision :: x()
    double precision, dimension(23):: x, y, z, w
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
    eps = 1e-12
    maxits = 479760

    do i = 1, n
        w(i) = 0
        x(i) = 1
        y(i) = i
        z(i) = i*i
    end do

    ! Método de Newton
    ! call metodo_newton(quadratica, quadratica_grad, quadratica_hess, x, alpha, beta, gama, sigma, rho, eps, maxits)
    ! call metodo_newton(rosenbrock, rosenbrock_grad, rosenbrock_hess, x, alpha, beta, gama, sigma, rho, eps, maxits)
    ! call metodo_newton(styblinsky_tang, styblinsky_tang_grad, styblinsky_tang_hess, x, alpha, beta, gama, sigma, rho, eps, maxits)
    ! call metodo_newton(rastrigin, rastrigin_grad, rastrigin_hess, x, alpha, beta, gama, sigma, rho, eps, maxits)

end program newton
