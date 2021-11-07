#=
Programa utilizado para gráfia, este programa utiliza las funciones de 
de la clase McCall para esto desarrrollare los ejercicioss

=#

## Exercise 1
# Parte 1
infb = 0.5;
supb = 1.5;
Nb = 25;

# β = range(infb, supb,length=Nb)
β = range(infb, supb,length=Nb);
ν = similar(β);



for jt in eachindex(β)
#    ν[jt] = vfi!(McCall(b=β[jt]))
    local mc  = McCall(b = β[jt],Nw=5000);
    vfi!(mc,verbose =false);
    ν[jt] = mc.w_star;
end

p1a = plot( β, ν,kind="scatter",mode="markers",Layout(
    title="Salario de reserva respecto a la compensacion de desempleo",
    xaxis_title="Compensacion de desempleo",
    yaxis_title="Salario de reserva"))
p1a

savefig(p1a, "PS1-g1a.png")
#Parte 2
inf_beta = 0.9;
sup_beta =0.99;
Nb = 25;

# β = range(infb, supb,length=Nb)
Β = range(inf_beta, sup_beta,length=Nb);
Α = similar(Β);



for jt in eachindex(Β)
#    ν[jt] = vfi!(McCall(b=β[jt]))
    local mc  = McCall(β = Β[jt],Nw=5000);
    vfi!(mc,verbose =false);
    Α[jt] = mc.w_star;
end

p1b = plot( Β, Α,kind="scatter",mode="markers",Layout(
    title="Salario de reserva respecto a la impaciencia del ejemplo",
    xaxis_title="Impaciencia",
    yaxis_title="Salario de reserva"))
    p1b

savefig(p1b, "PS1-g1b.png")

# p3 = scatter( β, ν)

## Exercise 2

mc =McCall();
vfi!(mc);
k =1000;
T = zeros(k)

T = dist_T(mc,k);

p2 = plot([histogram(x = T, opacity=0.75)],Layout(
    title="Distribucion de aceptar de oferta de salarios",
    xaxis_title="Periodos hasta que se acepto el salario",
    yaxis_title="Veces que se ofrecio el salario"))

savefig(p2, "PS1-g2.png")

#Estadisitcas descriptivas (Media, Mínimo, Máximo, Cunatiles =0.25,0.5,0.75)

print("\n Media :  \n")
print(mean(T))
print("\n Valor Mínimo : \n")
print(quantile(T,0))
print("\n Valor Máximo :  \n")
print(quantile(T,1))
print("\n Cuartil Inferior :  \n ")
print(quantile(T,0.25))
print("\n Mediana : \n")
print(quantile(T,0.5))
print("\n Cuartil superior :  \n")
print(quantile(T,0.75))


## Exercise 3

inf_beta = 0.9;
sup_beta = 0.99;
Nb = 100;

# β = range(infb, supb,length=Nb)
C = range(inf_beta, sup_beta,length=Nb);
ν = similar(C);



for jt in eachindex(C)
#    ν[jt] = vfi!(McCall(b=β[jt]))
    local mc  = McCall(β = C[jt],Nw=1000);
    vfi!(mc,verbose =false);
    ν[jt] = mean(dist_T(mc,10000));
end

# p3 = scatter( β, ν)


p3 = plot( C, ν,kind="scatter",mode="markers",Layout(
    title="Tiempo promedio para aceptar un salario dado un factor de descuento",
    xaxis_title="Factor de descuento",
    yaxis_title="Tiempo promedio"));

savefig(p3, "PS1-g3.png")

## Exercise 4

inf_theta = 0.1;
sup_theta = 10;
Nb = 100;

# β = range(infb, supb,length=Nb)
E5 = range(inf_theta, sup_theta,length=Nb);
ν = similar(E5);
period = similar(E5);


for jt in eachindex(E5)
    local mc3  = McCall(θ = E5[jt],Nw=5000);
    vfi!(mc3,verbose =false);
    ν[jt] = mc3.w_star;
    period[jt] = mean(dist_T(mc3,10000));
end


# p3 = scatter( β, ν)




p5 = plot(E5,period,kind="scatter",mode="markers",Layout(
    title="Tiempo promedio para aceptar un salario dado un factor de pesimismo",
    xaxis_title="Factor de pesimismo",
    yaxis_title="Tiempo promedio"));
p4 = plot( E5, ν,kind="scatter",mode="markers",Layout(
    title="Salario de reserva dado un factor de pesimismo ",
    xaxis_title="Factor de pesimismo",
    yaxis_title="Salario de reserva"));

savefig(p4, "PS1-g4.png")
savefig(p5, "PS1-g5.png")

p5
p3