program dfp_test
    use dfp_method
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    ! double precision :: x()
    double precision, dimension(24):: x, y, z, w
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
    eps = 1e-17
    maxits = 479760

    ! Definição do vetor inicial x0
    do i = 1, n
        x(i) = 1
        y(i) = 0
        z(i) = i
        w(i) = i**2
    end do

    ! Método DFP
    ! call metodo_dfp(quadratica, quadratica_grad, x, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_dfp(rosenbrock, rosenbrock_grad, y, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_dfp(styblinsky_tang, styblinsky_tang_grad, x, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_dfp(rastrigin, rastrigin_grad, x, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
end program dfp_test
