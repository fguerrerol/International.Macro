function adj_rand(N,z)
adj_1  = zeros(N,N);
nlinks = N*z/2;

while nlinks>0
   x= floor.(rand(1,2).*N) .+ 1;
   x = Int.(x)
   while adj_1[x[1],x[2]] == 1 || x[1] == x[2]
       x= floor.(rand(1,2).*N) .+ 1;
       x = Int.(x)

   end
   
   adj_1[x[1],x[2]]=1;
   adj_1[x[2],x[1]]=1;
   
   nlinks = nlinks - 1;
end
return adj_1

end


function Util(U,adj_copy)
    AgUtil = 0;
    healthy=0;
    unhealthy=0;
    for i =1:length(U)
        if U[i] == 0
           healthy =sum(adj_copy[i,:])*1.5;
        else
           unhealthy = sum(adj_copy[i,:])*0.5;
        end
    AgUtil =healthy - unhealthy + AgUtil;
    end

    return AgUtil/2;

end




function Decission(U,p,r,z,N,beta)
   # global adj,Restricted,beta,links_0,cap_testeo


    #print(Restricted)
    pre_result = zeros(1,length(U));
    U_noised = zeros(1,length(U));
    U_noised .= scramble(U);


    sel =  rand(1:N,1,trunc(Int64,floor(N*cap_testeo)));
    adj_temp = zeros(N,N)
    links_no_corta = flinks(adj);
    contagios_no_corta = (sum(U.*-(Restricted .-1)))*(1+p*z-r);
    Loss_no_corta = Loss(links_0, links_no_corta, beta, contagios_no_corta);
    for i = 1:length(sel)
        if Restricted[sel[i]] != 1
            if U_noised[sel[i]] != 0
                adj_temp .= adj;
                adj_temp[sel[i],:].=0;
                adj_temp[:,sel[i]].=0;
                links_corta = flinks(adj_temp);
                contagios_corta = (sum(U.*-(Restricted .-1)))*(1+p*z-r)-1-p*sum(adj[sel[i],:]);
                Loss_corta = Loss(links_0, links_corta, beta, contagios_corta);
                pre_result[sel[i]] = Loss_no_corta - Loss_corta;
            end
        end
    end
    return pre_result

end


function Loss(links_0,links,beta,contagios)    
    L = beta*((links_0 - links)^2) + (1-beta)*(contagios ^ 2);
    return L
end




function  sis_net_limit(N ::Int64, z::Number, nsteps::Int64, p0, r0, beta_i,test_cap)::Matrix{Float64}

    global adj,r,Restricted,links_0,adj_0,hosp,cap_testeo
    beta = beta_i;

    
    cap_testeo = test_cap;
    
    
    
    adj = adj_rand(N,z);
    adj_0 = copy(adj);
    flip    = rand(1:N,1,1)[1];
    U       = zeros(1,N);
    Restricted = zeros(1,N);

    
    U[flip] = 1;
    inf     = 1;
    rec     = 0;
    flag =0;
    tseries = zeros(nsteps,5);
    step_l    = 1;

    intervencion = trunc(Int64,floor(0.2 * nsteps));
 
    if step_l == 1
      links_0=flinks(adj_0);
    end
    while step_l <= nsteps
        if step_l == intervencion
            flag=1;
        end  
        U = infect_async(U,flag,N,p0,r0,z,beta);
        tseries[step_l,:] = [step_l,sum(U)/N,sum(Restricted)/N,sum(U.* - (Restricted .-1))/N,flinks(adj)/links_0];

        step_l += 1;
    end
 
    return tseries
end
 








function flinks(matrix)

    links = 0;
    
    for i =1:size(matrix)[1]
        links =sum(matrix[i,:]) + links;
    end
    links = links/2;
    return links
end








function scramble(U)
    error = 1;
    U_noised = zeros(1,length(U));
    
    for i= 1:length(U)
        if U[i]==1
            U_noised[i] = rand(1)[1] < error;
        else
            U_noised[i]= rand(1)[1] > error;
        end
    end
    return U_noised
end







function recollect(pre_result,total,U)
    #global Restricted, adj
    if sum(U) != 0        
        idx = partialsortperm(pre_result[1,:], 1:total, rev=true);
        for j =1:length(idx)
            if pre_result[idx[j]]>0
                adj[:,idx[j]] .= 0;
                adj[idx[j],:] .= 0;
                Restricted[idx[j]] = 1;
            end
        end 
    end
end



function Reconnect(U)
    #global adj,Restricted,adj_0
    
    for i = 1:length(U)
        if  sum(Restricted) <= 0.0000001* length(Restricted)
            return
        end
        if ((U[i] == 0) && (Restricted[i]==1))
            adj[:,i] = adj_0[:,i].* -(Restricted .- 1)' ;
            adj[i,:] = adj_0[i,:].* -(Restricted .- 1)' ;
            Restricted[i]=0;
            
        end
    end
end








function infect_async(U,flag,N,p,r,z,beta)
    #global adj,Restricted
    
    sel =  rand(1:N,1,N);
    
    if flag ===1 
        Reconnect(U);
        result = Decission(U,p,r,z,N,beta);   
        recollect(result,trunc(Int64,N/10),U);    
    end
 
    for i = 1:length(sel)
       x_0   = sel[i];
       
       if U[x_0] == 0     
           nei = findall(x->x==1,adj[x_0,:]);    
           inf = findall(x->x==1,U[nei]);     
           pp = rand(1,length(inf)) .< p; 
           if sum(pp) > 0
             U[x_0] = 1;
           end
       else
           if rand(1)[1] < r      
               U[x_0] = 0;
           end
       end
       
    end
    return U   
 end