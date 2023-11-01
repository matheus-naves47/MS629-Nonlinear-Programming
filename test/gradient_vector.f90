program gradient_vector
    !
    use gradient_method
    use functions
implicit none

double precision,dimension(5) :: z
integer,parameter::ndim=size(z)
z=(/1.,2.,3.,4.,5./)


print *, "Teste Derivada Numérica (Vetor Gradiente)"


print *, z
print *, grad(quadratica, z)

end program gradient_vector
