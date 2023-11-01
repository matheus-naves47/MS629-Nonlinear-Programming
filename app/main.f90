program OTISER
    ! Programa OTImização SEm Restrições
    use gradient_method ! Módulo do método do Gradiente
    ! use newtons_method ! Módulo do método de Newton Globalizado
    ! use cp1_method ! Módulo do método de Correção de Posto 1
    ! use dfp_method ! Módulo do método de Davidon, Fletcher e Powell
    use functions ! Módulo contendo as funções de teste

    implicit none

    ! Declaração de Parâmetros
    ! double precision :: x()
    double precision, dimension(6) :: x, y, z, w
    integer, parameter::n = size(x)
    double precision :: alpha ! Taxa de decréscimo suficiente na Regra de Armijo
    ! double precision :: beta ! Constante de proporcionalidade para a direção
    ! double precision :: gamma ! Constante relativa ao ângulo entre a direção e o gradiente
    double precision :: sigma ! Fator de diminuição na busca unidimensional
    ! double precision :: rho ! Incremento inicial na globalização do método de Newton
    double precision :: eps ! precisão para a norma do gradiente (critério de parada).
    integer :: i, maxits ! Contador e número máximo de iterações

    alpha = 1e-4
    sigma = 0.5
    eps = 1e-4
    maxits = 0

    ! Definição do vetor inicial x0
    do i = 1, n
        x(i) = 1
        y(i) = 0
        z(i) = 2*i + 1
        w(i) = 2*i
    end do

    ! Método do Gradiente
    call metodo_gradiente(quadratica, x, alpha, sigma, eps, maxits)
    call metodo_gradiente(rosenbrock, y, alpha, sigma, eps, maxits)
    call metodo_gradiente(styblinsky_tang, z, alpha, sigma, eps, maxits)
    call metodo_gradiente(rastrigin, w, alpha, sigma, eps, maxits)

end program OTISER
