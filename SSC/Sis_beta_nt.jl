include("Funciones.jl")
print("Inicio de programa")
using Optim, Statistics,Interpolations, Printf,
      LinearAlgebra, PlotlyJS, ColorSchemes, Distributions,Random
using MAT
N      = 300;
z      = 5.7;
r      = 0.6;
nsteps = 10000;
nave   =10;


#resultado = zeros(nsteps,5);
p = [0.7];
#p = collect(0.3:0.1:0.7);
#p = 0.5;
betas = [0 0.1 0.2 0.25 0.26 0.27 0.28 0.29 0.31 0.32 0.33 0.34 0.35 0.4 0.5 0.6 0.7 0.8 0.9 1];#betas = collect(0:0.2:1);
#betas = [0,0.5,1]
test_cap = [1.0];
#test_cap = [0.5];
#test_cap =0.05:0.05:1;
#test_cap= collect(0.2:0.2:1);
#test_cap =[0.05,0.1,0.25,0.4,0.65,0.7,0.9,1];

r=0.6;


p=0.7

cap_test=1

output = "Resultados_Beta"
function Simulacion_beta(output,N,z,nsteps,prob,r,betas,cap_test)
    economia_beta= zeros(nsteps, length(betas));
    contagios_beta = zeros(nsteps, length(betas));
    restringidos_beta = zeros(nsteps, length(betas));


    suma_economia = zeros(nsteps,nave);
    suma_contagios = zeros(nsteps,nave);
    suma_restringidos = zeros(nsteps,nave);
    for b = 1:length(betas)
        beta_i =betas[b]
        nave = 10;
        for i =1:nave
            #print("Iteracion   = $(i)")
            resultado = sis_net_limit(N, z, nsteps, prob, r, beta_i,cap_test,i);
            
            suma_economia[:,i] = resultado[:,5];
            suma_contagios[:,i] = resultado[:,2];
            suma_restringidos[:,i] = resultado[:,3];
            historia = resultado[:,1];
        end
        economia_beta[:,b] = mean(suma_economia, dims=2);
        contagios_beta[:, b] = mean(suma_contagios, dims=2);
        restringidos_beta[:, b] = mean(suma_restringidos, dims=2);
    end
    
    A =economia_beta;
    B =contagios_beta;
    C = restringidos_beta;
    
    matwrite(output, Dict(
                   "Economia" => A,
                   "Contagios" =>B ,
                   "Restrictos"=> C,
                    "Betas" => betas,
                    "Prob" => prob
                    ); compress = true)
end
