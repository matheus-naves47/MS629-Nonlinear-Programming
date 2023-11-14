program cp1_test
    use cp1_method
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    ! double precision :: x()
    double precision, dimension(24):: x,y,z,w
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

    ! Método CP1
    ! call metodo_cp1(quadratica, quadratica_grad, w, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_cp1(rosenbrock, rosenbrock_grad, w, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_cp1(styblinsky_tang, styblinsky_tang_grad, w, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""
    ! call metodo_cp1(rastrigin, rastrigin_grad, w, alpha, beta, gama, sigma, eps, maxits)
    ! print *, ""

end program cp1_test
