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
return adj_1,nlinks

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




function Decission(adj,U,Restricted,p,r,z,N,beta,cap_t,links_0)
    pre_result = zeros(1,length(U));
    U_noised = zeros(1,length(U));
    U_noised .= scramble(U);


    sel =  rand(1:N,1,trunc(Int64,floor(N*cap_t)));
    adj_temp = zeros(N,N);
    adj_temp .= adj;
    links_no_corta = flinks(adj);
    contagios_no_corta = (sum(U.*-(Restricted .-1)))*(1+p*z-r);
    Loss_no_corta = Loss(links_0, links_no_corta, beta, contagios_no_corta);
    for i = 1:length(sel)
        if Restricted[sel[i]] != 1
            if U_noised[sel[i]] != 0
                adj_temp[sel[i],:].=0;
                adj_temp[:,sel[i]].=0;
                links_corta = flinks(adj_temp);
                contagios_corta = (sum(U.*-(Restricted .-1)))*(1+p*z-r)-1-p*sum(adj[sel[i],:]);
                Loss_corta = Loss(links_0, links_corta, beta, contagios_corta);
                pre_result[sel[i]] = Loss_no_corta - Loss_corta;
                adj_temp[sel[i],:].=adj[sel[i],:];
                adj_temp[:,sel[i]].=adj[:,sel[i]];
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

    beta = beta_i;

    
    cap_t = test_cap;
    
    adj_0 = zeros(N,N)
    adj_m = zeros(N,N)
    adj_m,links_0 = adj_rand(N,z);
    adj_0 .= adj_m;
    flip    = rand(1:N,1,1)[1];
    U       = zeros(1,N);
    Restricted = copy(zeros(1,N));

    
    U[flip] = 1;

    tseries = zeros(nsteps,5);
    step_l    = 1;
    intervencion = 0;
    intervencion = trunc(Int64,floor(0.2 * nsteps));
    links_0 = N*z/2;
    local step_t = nsteps
    flag::Bool=false;
    for step_l in 1:step_t
        flag = step_l < 0.2*step_t ? false : true;
        
        tseries[step_l,:] = [step_l,sum(U)/N,sum(Restricted)/N,sum(U.* - (Restricted .-1))/N,flinks(adj_m)/links_0];

        adj_m,U,Restricted = infect_async(adj_m,U,Restricted,flag,N,p0,r0,z,beta,cap_t,links_0,adj_0);

    end
    

    return tseries
end
 








function flinks(adj)
    links = 0;
    longitud  =size(adj)[1]
    for i =1:longitud
        links =sum(adj[i,:]) + links;
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







function recollect(adj,Restricted,pre_result,total,U)
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
    return adj,Restricted
end



function Reconnect(adj,U,Restricted,adj_0)
    #global adj,Restricted,adj_0
    
    for i = 1:length(U)
        if  sum(Restricted) <= 0.0000001* length(Restricted)
            break
        end
        if ((U[i] == 0) && (Restricted[i]==1))
            adj[:,i] = adj_0[:,i].* -(Restricted .- 1)' ;
            adj[i,:] = adj_0[i,:].* -(Restricted .- 1)' ;
            Restricted[i]=0;
            
        end
    end
    return adj,Restricted
end








function infect_async(adj,U,Restricted,flag_0,N,p,r,z,beta,cap_t,links_0,adj_0)
    #global adj,Restricted
    
    sel =  rand(1:N,1,N);
    
    if flag_0 
        adj,Restricted = Reconnect(adj,U,Restricted,adj_0);
        result = Decission(adj,U,Restricted,p,r,z,N,beta,cap_t,links_0);   
        adj,Restricted =recollect(adj,Restricted,result,trunc(Int64,N/10),U);    
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
    return adj,U,Restricted   
 end