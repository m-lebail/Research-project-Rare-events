library(plotrix)
library(mvtnorm)
library(latex2exp) # permet d'ajouter des equations latex aux figures

library(ggplot2)


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

#N_testes = c(10,100,500,1000,2000,2500)
N_testes = c(3000)
#iters = c(5,10,20,30,40)
#epsilon = 0
T_iterations = 20
#N = 2000


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
      if(phi(X_cand,vect) > L_k){
        X_curr = X_cand
        accept = accept + 1
      }
    }
    
    chain[,i] = X_curr
  }
  
  
  acceptation = accept / (T_iterations-1)
  
  return(list(chain=t(chain),taux_accept=acceptation))
  
}



relative_error_plus_duration = sapply(N_testes, function(N){

     start_time <- proc.time()[[3]]
     
  #relative_error = sapply(1:essais, function(b){
    probas = apply(omega,2,function(vect){
      
      #vect = omega[,200]
      
      X = matrix(rnorm(d*N,0,1),nrow=d,ncol=N)
      
      p_0 = 0.25
      
      k = 1
      
      phi_X = phi(X,vect)
      
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
        
        phi_new = phi(new_vectors,vect)
        
        X[,phi_index_sort[2,num_col_remplacement]] = new_vectors
        
        vect_phi_insere = phi_index_sort[,num_col]
        if(length(num_col) == 1){
          vect_phi_insere = matrix(vect_phi_insere,nrow=d)
        }
    
        
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
      
      #print("sortie de la boucle for !")
      
      ##Nombre de particules au-dessus de tau
      N_L = sum(phi_index_sort[1,] > tau)
      
      L_k = tau
      
      ##numéro des vecteurs qui sont tels qu'évalués par phi, ils sont supérieurs à L_k
      num_col = (1:N)[phi_index_sort[1,] > L_k]
      
      ##numéro de colonnes des vecteurs à remplacer
      num_col_remplacement = (1:N)[-num_col]
      
      #sum(phi_index_sort[1,] == phi(X[,phi_index_sort[2,]]) )
      
      X_valides = X[,phi_index_sort[2,num_col]]  
      
      length_non_accepte = num_col[1] - 1
      
      if(length_non_accepte > 0){
        
        R = sample(num_col[1]:N,length_non_accepte,replace = T)
        
        chaines = lapply(R, function(r){
          X_0 = X[,phi_index_sort[2,r]]
          mec = MetropolisH(T_iterations,step,X_0,vect,L_k)
          mec$chain
        })
        
        #chaines = res[1,]
        
        new_vectors = sapply(1:length_non_accepte,function(i){
          t(chaines[[i]][T_iterations,])
        })
        
        X[,phi_index_sort[2,num_col_remplacement]] = new_vectors
      }
      
      p_hat = (N_L/N) * (p_0)^(k-1)  
      
      list(p_chap=p_hat,X=X) 
      
    })
    
    p_j = sapply(1:J,function(i){
      probas[[i]]$p_chap
    })
    
    p_j/pnorm(-tau)
    
    Samples_conditionnaly = lapply(1:J,function(i){
      probas[[i]]$X
    })
    
    #View(t(Samples_conditionnaly[[1]])%*%omega )
    
    sum_pj = sum(p_j)
    
    prob_vect = p_j / sum_pj
    
    #r = (N*J)*1/2*T_iterations/100
    
    r = N*J*(1/4)
    
    j = sample(1:J,r,replace=TRUE,prob=prob_vect)
    
    Xi = sapply(j,function(i){
      choix = sample(1:N,1)
      Samples_conditionnaly[[i]][,choix]
    })
    
     
    #Xi_bis = sapply(1:J, function(i){
    #  choix = sample(1:N,prob_vect[i]*r,replace=F)
    #  Samples_conditionnaly[[i]][,choix]
    #})
    
    #Xi_matrix_no_doublon = do.call("cbind",Xi_bis)
    
    matrix_wtx = t(Xi)%*%omega
    s = apply(matrix_wtx,1,function(x){
      sum(x > tau)
    })
    
    res_union = mean( sum_pj / s )
    
    up_bound = exp(-(tau^2)/2)
    
    re = abs(up_bound - res_union)/up_bound
  

    end_time <- proc.time()[[3]]
    
    duration = end_time - start_time
    
    matrix(c(re,duration),nrow=2)
      
})

relative_error_plus_duration

#esperance_re = apply(relative_error_mean,2,mean)
#esperance_re
#sd_re = apply(relative_error_mean,2,sd)
#sd_re




##N=500

#T_iterations = c(5,20,50,100)

iters = c(5,10,20,30,40)
relative_error_iters = c(0.020165371,0.008232904,0.001161281,0.001080973,0.001052272)

df_N_2000 = data.frame(iter=iters,re=relative_error_iters)

plt_N_2000 = ggplot(df_N_2000) +
  aes(x = log10(iter)) +
  geom_point(mapping = aes(y=log10(re)),color = 'red',size=3) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\log_{10}(Number of iterations)$'), 
       y = TeX('$\\log_{10}(relative error)$'),
       title = "Evolution of the relative error as a function of the number of 
       Kernel iterations at each step for N = 2000")

plt_N_2000


N_s  = c(10,100,500)#1000,2000,3000)
relative_error_Ns = c(0.1995129,0.1687747,0.01066514)#,0.00955331,0.004280787,0.002069835)

log10_df_iter_20 = data.frame(log10_N_s=log10(N_s),log10_re=log10(relative_error_Ns))

log10_plt_iter_20 = ggplot(log10_df_iter_20) +
  aes(x = log10_N_s) +
  geom_point(mapping = aes(y=log10_re),color = 'red',size=3) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\log_{10}$(Number of iterations)'), 
       y = TeX('$\\log_{10}$(relative error)'),
       title =TeX( "Evolution of the relative error as a function of the number of particles used for each probability $P_j$"))

log10_plt_iter_20

#Comparaison

Ns_intuit = c(100, 1000 , 5000 , 10000 , 50000 , 100000, 200000)

relative_error_Ns_intuit = c(0.604592,0.1936599 , 0.0802948 , 0.02928709  ,0.03046118, 0.02545173 , 0.0070901)

duration_intuit = c(0.890000 ,4.3200000, 33.4500000 ,83.83000000 ,1013.3 , 2827.06 ,11116.11)

duration_ALOE = c(10.0900000 , 95.9700000  , 516.99000000  ,1060.45  ,  2176.19  ,3348.32)


df_ALOE = data.frame(re_ALOE=relative_error_Ns,duration_ALOE=duration_ALOE)
df_intuit = data.frame(re_intuit=relative_error_Ns_intuit,duration_intuit=duration_intuit)

ggplot(df_ALOE) + aes(x = log(duration_ALOE),y=log(re_ALOE)) + geom_point(color = 'red',size=3) +
  geom_smooth(method = "lm", se = FALSE,col="red") +
  theme_bw() + geom_point(data=df_intuit,mapping = aes(x=log(duration_intuit),y=log(re_intuit)),color = 'blue',size=3) +
geom_smooth(data=df_intuit,mapping = aes(x=log(duration_intuit),y=log(re_intuit)),method = "lm", se = FALSE,col="blue") +
labs(x = TeX('$\\log_{10}$(computation time)'), 
     y = TeX('$\\log_{10}$(relative error)'),
     title =TeX( "Evolution in logarithmic scale of the relative error as a function of the computational time"))


#T_iterations
relative_error = relative_error_plus_duration[1,]
durations = relative_error_plus_duration[2,]


df_T_iter_10 = data.frame(N=N_testes,re=relative_error)


plt_T_iter_20 = ggplot(df_T_iter_20) +
  aes(x = N) +
  geom_point(mapping = aes(y=re),color = 'red',size=3) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = 'Number N of particles', 
       y = 'relative error',
       title = "Evolution of the relative error as a function of the number of particles N 
       for Kernel iterations equal to 20 at each step")

plt_T_iter_20
