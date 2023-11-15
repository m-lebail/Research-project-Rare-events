library(plotrix)
library(mvtnorm)
library(latex2exp) # permet d'ajouter des equations latex aux figures


rm(list=objects())
graphics.off()


d = 2

tau = 6

J = 360

essais = 10

#u_comp = rnorm(d*J)

#matrix_u = matrix(u_comp,nrow=d,ncol=J)
#norme = apply(matrix_u,2,function(x){
#  sqrt(sum(x^2))
#})
##to multiply each column by a value, we use diag matix with the vector and matrix multiplication
#omega = matrix_u%*%diag(1/norme)

#apply(omega,2,function(x){
#  sqrt(sum(x^2))
#})

sin_w = sin(2*pi*(1:J)/J)
cos_w = cos(2*pi*(1:J)/J)
omega = rbind(sin_w,cos_w)

phi = function(x,vect)
{
  t(vect)%*%x
}

max_phi = function(x)
{
  #values = apply(omega,2,phi,x=x)
  #max(values)
  max(phi(x,omega))
}



#N_testes = c(10,100,500,1000,2000)
#iters = c(5,20,50,100)
#epsilon = 0
T_iterations = 20
N_total = c(200000)


MetropolisH = function(T_iterations,step,X0,vect,L_k){
  
  chain = matrix(0,nrow=d,ncol=T_iterations)
  
  X_curr = X0
  
  accept = 0
  
  chain[,1] = X_curr
  
  for(i in 2:T_iterations){
    X_cand = (X_curr + step*matrix(rnorm(d,0,1),nrow=d,ncol=1)) #/ sqrt(1+step^2)
    
    log_r = -(1/2) * (sum(X_cand^2) -  sum(X_curr^2) ) 
    
    u = runif(1)
    
    if(log(u) < log_r){
      if(max_phi(X_cand) > L_k){
        X_curr = X_cand
        accept = accept + 1
      }
    }
    
    chain[,i] = X_curr
  }
  
  
  acceptation = accept / (T_iterations-1)
  
  return(list(chain=t(chain),taux_accept=acceptation))
  
}



re_plus_duration = sapply(N_total, function(N){
 
    start_time <- proc.time()[[3]]  
  
    #vect = omega[,1]  
    
    X = matrix(rnorm(d*N,0,1),nrow=d,ncol=N)
    
    p_0 = 0.25
    
    k = 1
    
    phi_X = apply(X,2,max_phi)
    
    phi_index = rbind(phi_X,1:N)
    
    phi_index_sort = phi_index[,order(phi_index[1,],decreasing = F)]
    
    L_k = quantile(phi_index_sort[1,],1-p_0)
    
    
    step = 0.35
    
    
    
    while(L_k < tau){
      
      ##numéro des vecteurs qui sont tels qu'évalués par phi, ils sont supérieurs à L_k
      num_col = (1:N)[phi_index_sort[1,] > L_k]
      
      ##numéro de colonnes des vecteurs à remplacer
      num_col_remplacement = (1:N)[-num_col]
      
      #sum(phi_index_sort[1,] == phi(X[,phi_index_sort[2,]]) )
      
      X_valides = X[,phi_index_sort[2,num_col]]  
      
      length_non_accepte = num_col[1] - 1
      
      R = sample(num_col[1]:N,length_non_accepte,replace = T)
      
      chaines = lapply(R, function(r){
        X_0 = X[,phi_index_sort[2,r]]
        mec = MetropolisH(T_iterations,step,X_0,vect,L_k)
        mec$chain
      })
      
      #taux_acceptations = res[2,]
      #unlist(taux_acceptations)
      
      #chaines = res[1,]
      
      new_vectors = sapply(1:length_non_accepte,function(i){
        t(chaines[[i]][T_iterations,])
      })
      
      phi_new = apply(new_vectors,2,max_phi)
      
      X[,phi_index_sort[2,num_col_remplacement]] = new_vectors
      
      vect_phi_insere = phi_index_sort[,num_col]
      
      l = N-length_non_accepte
      
      
      
      for(i in 1:length_non_accepte){
        
        spot = findInterval(phi_new[i],vect_phi_insere[1,])
        if(1 <= spot && spot < l){
          vect_phi_insere = cbind(  vect_phi_insere[,1:(spot)] ,  rbind(phi_new[i],phi_index_sort[2,num_col_remplacement[i]] ) , vect_phi_insere[,(spot+1):l] )
        } else if(spot == l) {
          vect_phi_insere = cbind(  vect_phi_insere[,1:l] ,  rbind(phi_new[i],phi_index_sort[2,num_col_remplacement[i]] ))
        } else{ #spot == 0
          vect_phi_insere = cbind(rbind(phi_new[i],phi_index_sort[2,num_col_remplacement[i]] ),vect_phi_insere[,1:l]) 
        }
        
        l = l + 1
      }
      
      phi_index_sort = vect_phi_insere
      
      #sum(phi_index_sort[1,] == phi(X[,phi_index_sort[2,]]) )
      
      L_k = quantile(phi_index_sort[1,],1-p_0)
      
      k = k + 1
    }
    
    ##Nombre de particules au-dessus de tau
    N_L = sum(phi_index_sort[1,] > tau)
    
    #L_k = tau
    
    ##numéro des vecteurs qui sont tels qu'évalués par phi, ils sont supérieurs à L_k
    #num_col = (1:N)[phi_index_sort[1,] > L_k]
    
    ##numéro de colonnes des vecteurs à remplacer
    #num_col_remplacement = (1:N)[-num_col]
    
    #sum(phi_index_sort[1,] == phi(X[,phi_index_sort[2,]]) )
    
    #X_valides = X[,phi_index_sort[2,num_col]]  
    
    #length_non_accepte = num_col[1] - 1
    
    #if(length_non_accepte > 0){
      
      #R = sample(num_col[1]:N,length_non_accepte,replace = T)
      
      #chaines = lapply(R, function(r){
        #X_0 = X[,phi_index_sort[2,r]]
        #mec = MetropolisH(T_iterations,step,X_0,vect,L_k)
        #mec$chain
      #})
      
      #chaines = res[1,]
      
      #new_vectors = sapply(1:length_non_accepte,function(i){
      #  t(chaines[[i]][T_iterations,])
      #})
      
      #X[,phi_index_sort[2,num_col_remplacement]] = new_vectors
    #}
    
    
    p_hat = (N_L/N) * (p_0)^(k-1)  
    
    up_bound = exp(-(tau^2)/2)
  
    re = abs(up_bound - p_hat)/up_bound
    
    end_time <- proc.time()[[3]]
    
    duration = end_time - start_time
    
    matrix(c(re,duration),nrow=2)
    
})

    
re_plus_duration  

 esperance_re = mean(relative_error)
esperance_re
sd_re = sd(relative_error)
sd_re


