library(plotrix)
library(mvtnorm)
library(latex2exp) # permet d'ajouter des equations latex aux figures


rm(list=objects())
graphics.off()


d = 2

tau = 6

J = 360

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


T_iterations = 20
N = 50000

start_time <- Sys.time()

m = 1

##We try first with a ditribution mu equals to standard gaussian distribution

X = matrix(rnorm(d*N,0,1),nrow=d,ncol=N)

phi_X = max_phi(X)

phi_index = rbind(phi_X,1:N)

phi_index_sort = phi_index[,order(phi_index[1,],decreasing = F)]

L_m = phi_index_sort[1,1]



#top = 267459
while(L_m < tau){
  #for(j in 1:top){  
  R = sample(2:N,1)
  
  X_m_plus_un = X[1:d,phi_index_sort[2,R]]
  
  #T_iterations = 30
  
  ##Metropolis-Hastings sampler with Gaussiand random walk with a supplementary condition
  #chain = matrix(0,nrow=d,ncol=T_iterations)
  step = 0.35
  #sigma = diag(d)
  
  X_curr = matrix(X_m_plus_un,nrow=d,ncol=1)
  #chain[,1] = X_curr 
  
  for(i in 2:T_iterations){
    X_cand = (X_curr + step*matrix(rnorm(d,0,1),nrow=d,ncol=1)) 
    
    log_r = -(1/2) * (sum(X_cand^2) -  sum(X_curr^2) )
    
    #r = exp(-(1/2) * (sum(X_cand^2) -  sum(X_curr^2) ))
    
    u = runif(1)
    
    if(log(u) < log_r){
      if(max_phi(X_cand) > L_m ){
        X_curr = X_cand
      }
    }
    
    #chain[,i] = X_curr
  }
  
  phi_new = max_phi(X_curr)
  
  X[,phi_index_sort[2,1]] = X_curr
  
  spot = findInterval(phi_new, phi_index_sort[1,2:N])
  
  if(1 <= spot && spot < N-1){
    phi_index_sort = cbind(  phi_index_sort[,2:(spot+1)] ,  rbind(phi_new,phi_index_sort[2,1]) , phi_index_sort[,(spot+2):N] )
  } else if(spot == N-1) {
    phi_index_sort = cbind(  phi_index_sort[,2:N] ,  rbind(phi_new,phi_index_sort[2,1]))
  } else{ #spot == 0
    phi_index_sort = cbind(rbind(phi_new,phi_index_sort[2,1]),phi_index_sort[,2:N]) 
  }
  
  L_m = phi_index_sort[1,1]
  
  m = m + 1
}



p_chap = (1 - (1/N))^(m-1)

up_bound = exp(-(tau^2)/2)

re = abs(up_bound - p_hat)/up_bound


stop_time <- Sys.time()
duration = stop_time - start_time
duration

re
