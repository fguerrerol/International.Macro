include("Funciones.jl")
using Optim, Statistics,Interpolations, Printf, LinearAlgebra, PlotlyJS, ColorSchemes, Distributions
using MAT
N      = 300;
z      = 5.7;
r      = 0.6;
nsteps = 10000;
nave   =10;

#resultado = zeros(nsteps,5);

p = collect(0.3:0.1:0.7);

betas = collect(0:0.2:1);
#test_cap =0.05:0.05:1;
test_cap= collect(0.3:0.1:1);
#test_cap =[0.05,0.1,0.25,0.4,0.65,0.7,0.9,1];

r=0.6;
economia_segun_beta = zeros(nsteps, length(test_cap),length(betas));
contagios_segun_beta= zeros(nsteps, length(test_cap),length(betas));
restringidos_segun_beta= zeros(nsteps, length(test_cap),length(betas));

econ_prob = zeros(nsteps, length(test_cap),length(betas),length(p));
contagios_prob = zeros(nsteps, length(test_cap),length(betas),length(p));

restringidos_prob = zeros(nsteps, length(test_cap),length(betas),length(p));


economia_test= zeros(nsteps, length(test_cap));
contagios_test = zeros(nsteps, length(test_cap));
restringidos_test = zeros(nsteps, length(test_cap));


suma_economia = zeros(nsteps,nave);
suma_contagios = zeros(nsteps,nave);
suma_restringidos = zeros(nsteps,nave);


for c = 1:length(p)
  prob = p[c];
  
  for k = 1:length(betas)
   
    beta_i = betas[k];



    for j = 1:length(test_cap)  
         
      
      Threads.@threads for i =1:nave
          #print("Iteracion   = $(i)")
          resultado = sis_net_limit(N, z, nsteps, prob, r, beta_i,test_cap[j]);
          
          suma_economia[:,i] = resultado[:,5];
          suma_contagios[:,i] = resultado[:,2];
          suma_restringidos[:,i] = resultado[:,3];
          historia = resultado[:,1];
      end
      #print(j)
      economia_test[:,j] = mean(suma_economia, dims=2);
      contagios_test[:, j] = mean(suma_contagios, dims=2);
      restringidos_test[:, j] = mean(suma_restringidos, dims=2);
    end
      economia_segun_beta[:, :,k] = economia_test;
      contagios_segun_beta[:,:, k] = contagios_test;
      restringidos_segun_beta[:,:, k] = restringidos_test;
      
  end
  econ_prob[:,:,:,c] =economia_segun_beta;
  contagios_prob[:,:,:,c] = contagios_segun_beta;
  restringidos_prob[:,:,:,c] =restringidos_segun_beta;

end

A =econ_prob;
B =contagios_prob;
C = restringidos_prob;

matwrite("ResultadosprobF.mat", Dict(
               "Economia" => economia_segun_beta,
               "Contagios" =>contagios_segun_beta ,
               "Restrictos"=> restringidos_segun_beta,
                "Betas" => betas,
                "probabilitis" => p,
                "test_cap" => test_cap
                ); compress = true)