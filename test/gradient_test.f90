program gradiente
    use gradient_method
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    ! double precision :: x()
    double precision, dimension(23):: x, y, z, w
    integer, parameter::n = size(x)
    double precision :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
    double precision :: sigma ! Fator de diminuição na busca unidimensional
    double precision :: eps ! precisão para a norma do gradiente (critério de parada).
    integer :: i, maxits ! Contador e número máximo de iterações

    alpha = 1e-4
    sigma = 0.5
    eps = 1e-12
    maxits = 479760

    ! Definição do vetor inicial x0
    do i = 1, n
        w(i) = 0
        x(i) = 1
        y(i) = i
        z(i) = i*i
    end do

    ! Método do Gradiente

    call metodo_gradiente(quadratica, quadratica_grad, x, alpha, sigma, eps, maxits)
    ! call metodo_gradiente(rosenbrock, rosenbrock_grad, x, alpha, sigma, eps, maxits)
    ! call metodo_gradiente(styblinsky_tang, styblinsky_tang_grad, x, alpha, sigma, eps, maxits)
    ! call metodo_gradiente(rastrigin, rastrigin_grad, x, alpha, sigma, eps, maxits)

end program gradiente
