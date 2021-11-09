include("Funciones.jl")
using Optim, Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using MAT
N      = 300;
z      = 5.7;
r      = 0.6;
nsteps = 500;
nave   =6;

resultado = zeros(nsteps,5);

p = [0.5];

betas = [1];
test_cap =[0.7,0.9,1]
#test_cap =[0.05,0.1,0.25,0.4,0.65,0.7,0.9,1];
prob = p[1];
r=0.6;
economia_segun_beta = zeros(nsteps, length(test_cap),length(betas));
contagios_segun_beta= zeros(nsteps, length(test_cap),length(betas));
restringidos_segun_beta= zeros(nsteps, length(test_cap),length(betas));


for k = 1:length(betas)
  economia_test= zeros(nsteps, length(test_cap));
  contagios_test = zeros(nsteps, length(test_cap));
  restringidos_test = zeros(nsteps, length(test_cap));
  beta_i = betas[k];



  Threads.@threads for j = 1:length(test_cap)  
    
    i=1;
    
    output = "Probabilidad de contagio  = $(prob) ,/n T.Recuperación r = $(r),
     Ponderador del gobierno Beta = $(beta_i), /n Cap.Testeo = $(test_cap[j])";
    print(output)
    
    
    suma_economia = zeros(nsteps);
    suma_contagios = zeros(nsteps);
    suma_restringidos = zeros(nsteps);
    while i<=nave
        #print("Iteracion   = $(i)")
        global resultado = sis_net_limit(N, z, nsteps, 0.5, r, beta_i,test_cap[j]);
        
        suma_economia = resultado[:,5] + suma_economia;
        suma_contagios = resultado[:,2] + suma_contagios;
        suma_restringidos = resultado[:,3] + suma_restringidos;
        historia = resultado[:,1];

        i=i+1;  
    end
    print(j)
    economia_test[:,j] = suma_economia./nave;
    contagios_test[:, j] = suma_contagios./nave;
    restringidos_test[:, j] = suma_restringidos./nave;
  end
    economia_segun_beta[:, :,k] = economia_test;
    contagios_segun_beta[:,:, k] = contagios_test;
    restringidos_segun_beta[:,:, k] = restringidos_test;
    
end

A =economia_segun_beta;
B =contagios_segun_beta;
C = restringidos_segun_beta;

matwrite("Resultados_julia2.mat", Dict(
               "Economia" => economia_segun_beta,
               "Contagios" =>contagios_segun_beta ,
               "Restrictos"=> restringidos_segun_beta,
                "Betas" => betas,
                "Cap_test" =>test_cap 

                   ); compress = true)