using Statistics

include("SGU.jl")
#### Modelo normal
A = SOEwr();
comp_eqm!(A)

### Modelo con el valor más bajo de w 

p3_1a = plot(contour(x =A.agrid,
                     y =A.zgrid,
                     z = A.Y,
                     colorscale="inferno",
                     contours_start = 0, 
                     contours_size = 0.05),
                     Layout(tite ="Curvas de nivel de producto (wbar =0)",
                     xaxis_title="Ahorro",
                     yaxis_title="Productividad"))

savefig(p3_1a, "ps3-p1_1A.png")


B = SOEwr(wbar = 0.0001);
comp_eqm!(B)

p3_1b = plot(contour(x =B.agrid,
                     y =B.zgrid,
                     z = B.Y,
                     colorscale="inferno",
                     contours_start = 0, 
                     contours_size = 0.05),
                     Layout(tite ="Curvas de nivel de prodcuto (wbar =0)",
                     xaxis_title="Ahorro",
                     yaxis_title="Productividad"))

savefig(p3_1b, "ps3-p1_1B.png")
#### Para simular utilizo la funcion simul y le hago 200 periodos de simulacion
path_a = simul(A,T=  250)

path_b = simul(B,T = 250)

#### Leugo de simular hallo la volatilidad relativa
#### usando la desviaciṕn estándar tanto de el comsumo y 
### el producto

##### Primera volatitlidaf

std(path_a[:C])/std(path_a[:Y])




#### Seungda volatildiad 

std(path_b[:C])/std(path_b[:Y])





#### Covariannza
#### Covarianza del primer elemento


cov(path_a[:Y],path_a[:CA])

### Covarianza del segundo

cov(path_b[:Y],path_b[:CA])



p1= plot([histogram(x =path_a[:A], opacity=0.75)], Layout(title="Distribución de la deuda caso default"))


p2= plot([histogram(x =path_b[:A], opacity=0.75)], Layout(title="Distribución de la deuda caso default"))




savefig(p1, "Distrib_ergo.png")
savefig(p2, "Dsitrib_ergoB.png")