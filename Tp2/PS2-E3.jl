using StatsBase
using LatexPrint
include("arellano.jl")
include("ifp.jl")
include("def_simul.jl")
include("cakeeating.jl")
include("itpcake.jl")
cd("/home/francisco/Escritorio/Open-macro/Francisco/")



## Primer Ejercicio de Problem Set

ce = CakeEating();
vfi_itp!(ce, verbose = true,tol=1e-10);
PS2_graph1 = plot(scatter(x = ce.kgrid, y = ce.gc), Layout(title = "Consumo en función del capital"))
savefig(PS2_graph1, "gráifco_1_PS2.png");
#Este gráfico es una recta que simboliza el consumo

ratio_ck = ce.gc ./ ce.kgrid  
PS2_graph2 = plot(scatter(x = ce.kgrid, y = ratio_ck), Layout(title = "Ratio c/k en función de k"))
savefig(PS2_graph2, "gráfico_2_PS2.png")

### Cuando el capital es menor a 0.10



ratio_kk = ce.gc ./ ce.kgrid
PS2_graph3 = plot(scatter(x = ce.gc, y = ratio_kk), Layout(title = "Ratio c/k en función de k"))
savefig(PS2_graph2, "gráfico_2_PS2.png")



#Mostrar la funcion de ahorro

p1 = scatter(x = ce.kgrid, y =ce.gk./ce.kgrid,name="Consumo")
p2 = scatter(x = ce.kgrid, y =ce.gc./ce.kgrid,name="Ahorro")

p=[p1 p2]


PS2_graph4 = plot([p1,p2], Layout(title = "Consumo y ahorro en función del capital",xaxis_title="Capital"))
savefig(PS2_graph4, "Grafico1.png");


#### Simualcion Torta #### 
function Torta(;T::Number=50,k_0::Number=1.0)
    C = zeros(T);
    K = zeros(T);

    if k_0 != 1.0
       local ce =CakeEating(k = k_0,Nk = 200);
    else
       local ce =CakeEating(Nk=200);
    end
    vfi_itp!(ce);
    k_0 = 0;
  

    itp_v = interpolate((ce.kgrid,), ce.gc, Gridded(Linear()));


    for i in 1:T 
        if i == 1
            
            k_0 = maximum(ce.kgrid)
            c = itp_v(maximum(ce.kgrid))
            k_1 = k_0  -c;
        else
            c =itp_v(k_0*(1+ce.r))
            k_1 = k_0*(1+ce.r) -c;
        end
        k_0= k_1;
        C[i] = c;
        K[i] = k_1;
    end
    plot_consumo=scatter(x=T,y=C,name="Consumo")
    plot_torta=scatter(x=T,y=K,name="Ahorro")
    A = plot([plot_consumo, plot_torta], Layout(title="Consumo y torta restante", xaxis_title="Tiempo"))
    return A
end

Torta_graph = Torta(T=400)
savefig(Torta_graph,"Evolution-torta.png")



#### Segundo Ejercico de Problem Set
c = IFP()
vfi!(c)
c.
PS2_p21 =   plot(contour(x=c.kgrid,y=c.ygrid,z= c.gc,colorscale="inferno"),
Layout(title="Consumo dado nivel de capital y estado de la economia", yaxis_title="Estado de la economía",xaxis_title="Nivel de capital"))


PS2_p2 = plot([scatter(x = c.kgrid, y = c.gc[:,jy],name= "y = $(c.ygrid[jy])") for jy in eachindex(c.ygrid)])
savefig(PS2_p21, "Diferentes niveles de ingreso.png")


mutable struct Elementos
	Con::Vector{Float64}
	Kap::Vector{Float64}	
	Yield::Vector{Float64}
end


function simul_ingresop1(;k_0 = 1,θ = 0, tol = 1e-5)
    if k_0 != 1.0
       global ingreso =IFP(k = k_0,Nk = 100,θ=θ);
     else
       global ingreso =IFP(Nk=100,θ=θ);
     end
     vfi!(ingreso,tol=tol);
end

function simul_ingresos(;T::Number=40)
    C = zeros(T);
    K = zeros(T);
    Y = zeros(T);
    global Time=1:1:T
    k_0 = 0

  
    knots = (ingreso.kgrid,ingreso.ygrid)
    itp_c = interpolate(knots, ingreso.gc, Gridded(Linear()));
    itp_k = interpolate(knots, ingreso.gk, Gridded(Linear()));
    y_index=0

    Ny = 1:1:25

    for i in Time 
        if i == 1
            y_index = 13;
            k_0 = maximum(ingreso.kgrid)                   
        else
            y_index = sample(Ny,Weights(ingreso.Py[:,y_index]))
        end
        c = itp_c(k_0,ingreso.ygrid[y_index])
        k_1 = itp_k(k_0,ingreso.ygrid[y_index]) 
        k_0= k_1;
        C[i] = c;
        K[i] = k_1/(1+ingreso.r);
        Y[i] = ingreso.ygrid[y_index];

    end
    plot_consumo=scatter(x=Time,y=C,name="Consumo")
    plot_ahorro=scatter(x=Time,y=K,name="Ahorro")
    plot_dotacion=scatter(x=Time,y=Y,name="Dotacion")

    Evolution = plot([plot_consumo, plot_ahorro,plot_dotacion], Layout(title="Consumo y Ahorro restante", xaxis_title="Tiempo"))
    savefig(Evolution,"Evolution.png")
    global Grafo = Elementos(C,K,Y);
    
    p1= plot([histogram(x = Grafo.Con./Grafo.Yield, opacity=0.75)], Layout(title="Distribución de c/y"))
    global p2 = plot([histogram(x = Grafo.Kap, opacity=0.75)], Layout(title="Distribución del capital"))
    return p1
end


simul_ingresop1(θ = 0)
simul_ingresos()
norob =Grafo
simul_ingresop1(θ = 10)
simul_ingresos()
rob10 = Grafo
simul_ingresop1(θ = 30)
simul_ingresos()
rob30 = Grafo
simul_ingresop1(θ = 100)
simul_ingresos()
rob100 = Grafo

no_rob=scatter(x=Time,y=norob.Con./norob.Yield,name="Dotacion")
rob10=scatter(x=Time,y=rob10.Con./rob10.Yield,name="Dotacion")
rob_30=scatter(x=Time,y=rob30.Con./rob30.Yield,name="Dotacion")
rob_100=scatter(x=Time,y=rob100.Con./rob100.Yield,name="Dotacion")


Propension = plot([no_rob,rob_10, rob_30,rob_100], Layout(title="Consumo y Ahorro restante", xaxis_title="Tiempo"))




#### Tercer Ejercicio de Problem Set

A = Default();
mpe!(A,tol=1e-5);

C = Default(Δ = 0.08);
mpe!(C,tol=1e-5);

D = Default(ℏ = 0.7);
mpe!(D,tol=1e-5);

B = Default(β = 0.98);
mpe!(B,tol=1e-5);


A_graph = plot([scatter(x = A.bgrid, y = A.prob[:,jy],name= "y = $(A.ygrid[jy])") for jy in eachindex(A.ygrid)],Layout(title = "Párametros normales"))
B_graph = plot([scatter(x = B.bgrid, y = B.prob[:,jy],name= "y = $(B.ygrid[jy])") for jy in eachindex(B.ygrid)],Layout(title = "Mayor factor de descuento"))


C_graph = plot([scatter(x = C.bgrid, y = C.prob[:,jy],name= "y = $(C.ygrid[jy])") for jy in eachindex(C.ygrid)],Layout(title ="Mayor coste de Default"))

D_graph = plot([scatter(x = D.bgrid, y = D.prob[:,jy],name= "y = $(D.ygrid[jy])") for jy in eachindex(D.ygrid)],Layout(title="Mayor haircut"))






savefig(B_graph, "Mayor factor de descuento.png")
savefig(A_graph, "Parametros normales.png")
savefig(C_graph, "Menor Coste de Default.png")
savefig(D_graph, "Mayor haircut.png")

### 3.2.1 Distribucion de la Deuda

Deuda_corta , Deuda_larga = samples_comp()

Combined_graph = plots_comp(Deuda_corta, Deuda_larga)


savefig(Combined_graph, "Distribucion deuda.png")


### 3.2.2 Correlaciones

var_consumption = var(Deuda_larga[:c])

var_product = var(Deuda_larga[:y])

sdivy = cor(Deuda_larga[:y],Deuda_larga[:spread])


Tabla = Array{Any}(nothing, 3, 2);

Tabla[1,1] = "Variance of comsuption";
Tabla[2,1] = "Variance of product";
Tabla[3,1] = "Correlation between spread and product";
Tabla[1,2] = var_consumption;
Tabla[2,2] = var_product;
Tabla[3,2] = sdivy;

tabular(Tabla, hlines=true);


# 3.2.3 Robustness
short, long = samples_comp()
rob_25_short, rob_25_long = samples_comp(ι = 25)
rob_100_short, rob_100_long = samples_comp(ι = 100)

r_dy_short_rob_25 = rob_25_short[:d]./rob_25_short[:y]
r_dy_long_rob_25 = rob_25_long[:d]./rob_25_long[:y]

r_dy_short = short[:d]./short[:y]
r_dy_long = long[:d]./long[:y]

r_dy_short_rob_100 = rob_100_short[:d]./rob_100_short[:y]
r_dy_long_rob_100 = rob_100_long[:d]./rob_100_long[:y]
### 3.3 Deuda contingente

#### Fluctuacion de ingresos ####
# Inicializamos el modelo104

NoCon = Default();
mpe!(NoCon, tol=1e-5)

Contin = Default(α = 0.4);
mpe!(Contin, tol =1e-5)

p4b =   plot(contour(x=NoCon.bgrid,y=NoCon.ygrid,z= (-1)*NoCon.v,colorscale="inferno"),
Layout(title="Funcion de valor dado ingreso y dada deuda", xaxis_title="Deuda",yaxis_title="Ingreso"))


savefig(p4b, "ps2-p4b.png")

p5b =   plot(contour(x=Contin.bgrid,y=Contin.ygrid,z= (-1)*Contin.v,colorscale="inferno"),
Layout(title="Funcion de valor dado ingreso y dada deuda", xaxis_title="Deuda",yaxis_title="Ingreso"))

savefig(p5b, "ps2-p5b.png")


p231c = contour(x=NoCon.bgrid,y=NoCon.ygrid,z= (-1)*NoCon.v,colorscale="inferno",name="Deuda No contingente")

p232c = contour(x=Contin.bgrid,y=Contin.ygrid,z= (-1)*Contin.v,colorscale="inferno",name="Deuda Contingente")

pfin = plot([p231c, p232c],Layout(barmode="overlay"))
