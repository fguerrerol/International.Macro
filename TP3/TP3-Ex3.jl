using JLD2
using Optim, Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using QuantEcon: tauchen

include("SGU_ps3_v3.jl")




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
