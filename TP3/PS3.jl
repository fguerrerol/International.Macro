using JLD2, Statistics
using Optim, Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using QuantEcon: tauchen


##### Parte 1 ########


include("SGU.jl")
include("SGU_cc.jl")


A = SOEwr();
comp_eqm!(A)

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



C = SOEwr_cc();
comp_eqm!(C)


save_object(C,"SOE_cc.jdl2");


C = load_object("SOE_cc.jdl2");


##### Grafico los activos
p1 =plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Y[1,:,:],colorscale = "inferno"),
Layout(title="Nivel de producto dado shock de producividad y shock de noticias bajo
        endeudamiento",
            xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p1,"Endeudados.png")

p2 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Y[18,:,:],colorscale="inferno"),
Layout(title="Nivel de producto dado shock de producividad y shock de noticias bajo
    activos casi nulos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p2,"Cero-deuda.png")

p3 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Y[40,:,:],colorscale="inferno"),
Layout(title="Nivel de producto dado shock de producividad y shock de noticias bajo
    el nivel más alto de activos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))


savefig(p3,"Mayor cantidad de activos.png")




##### Para el ahorro del siguiente periódo


##### Grafico los activos
p4 =plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Ap[1,:,:],colorscale = "inferno"),
Layout(title="Ahorro próximo shock de producividad y shock de noticias bajo
endeudamiento",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p4,"InversionA0.png")

p5 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Ap[20,:,:],colorscale="inferno"),
Layout(title="Ahorro próximo dado shock de producividad y shock de noticias bajo
activos casi nulos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p5,"InversionA18.png")

p6 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= C.Ap[40,:,:],colorscale="inferno"),
Layout(title="Ahorro próximo dado shock de producividad y shock de noticias bajo
el nivel más alto de activos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))


savefig(p6,"InversionA40.png")


consumo = zeros(40,15,15)
for i = 1:15
consumo .= C.gc[i,:,:,:] + consumo
end

consumo_f = consumo/15

#### Graficos de la cuenta corriente

neg_asset = (C.Y - consumo_f)[1,:,:]

zero_asset = (C.Y- consumo_f)[18,:,:]

pos_asset = (C.Y-consumo_f)[40,:,:]

##### Grafico los activos
p7 =plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= neg_asset,colorscale = "inferno"),
Layout(title="Cuenta corriente shock de producividad y shock de noticias bajo
endeudamiento",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p7,"CurrentA1.png")

p8 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= zero_asset,colorscale="inferno"),
Layout(title="Cuenta corriente dado shock de producividad y shock de noticias bajo
activos casi nulos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))

savefig(p8,"CurrentA18.png")

p9 = plot(contour(x =C.zgrid, 
y = C.ξgrid,
z= pos_asset,colorscale="inferno"),
Layout(title="Cuenta Corriente dado shock de producividad y shock de noticias bajo
el nivel más alto de activos",
xaxis_title="Shock de productividad",yaxis_title="Shock de noticias"))


savefig(p9,"CurrentA40.png")
