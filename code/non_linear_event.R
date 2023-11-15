library(plotrix)
library(mvtnorm)
library(latex2exp) # permet d'ajouter des equations latex aux figures
library(misc3d)
library(plotrix)

rm(list=objects())
graphics.off()

J = 2


sin_w = sin(2*pi*(1:J)/J)
cos_w = cos(2*pi*(1:J)/J)
omega = rbind(sin_w,cos_w)

#u = matrix(omega[,1],nrow=2)

phi = function(X,u)
{
  norme = sqrt(t(X)%*%X)
  return( abs(t(X)%*%u) / norme )
}

frontiere = function(x,y,u)
{
  X = rbind(x,y)
  apply(X,2,phi,u=u)
}

d = 2

#J = 2


x = seq(-tau-0.2,tau+0.2,length.out=100)
y = seq(-tau-0.2,tau+0.2,length.out=100)




tau = 0.95

par(pty="s") #square plot region

plot(0,0,xlim=c(-tau-0.5,tau+0.5),ylim=c(-tau-0.5,tau+0.5))

draw.circle(0,0,tau)

apply(omega,2, function(u){
  z = outer(x, y, frontiere,u=u)
  contour(x,y,z,levels=tau,add = T,col="green",lwd=2)
})



f = function(x,y,theta)
{
  x_bis = x*cos(theta) - y*sin(theta)
  y_bis = x*sin(theta) + y*cos(theta)
  (1/20)*(y_bis^2) + x_bis  
}



par(pty="s") #square plot region

plot(0,0,xlim=c(-4.5,4.5),ylim=c(-4.5,4.5))

draw.circle(0,0,1)

x = seq(-4,4,length.out=100)
y = seq(-4,4,length.out=100)

theta = 0
z = outer(x, y, f,theta=theta)
contour(x,y,z,levels=3,add = T,col="red",lwd=2)

theta = pi/2
z = outer(x, y, f,theta=theta)
contour(x,y,z,levels=3,add = T,col="red",lwd=2)


theta = pi
z = outer(x, y, f,theta=theta)
contour(x,y,z,levels=3,add = T,col="red",lwd=2)

theta = 3*pi/2
z = outer(x, y, f,theta=theta)
contour(x,y,z,levels=3,add = T,col="red",lwd=2)


phi = function(Xvect,theta)
{
  x = Xvect[1,]
  y = Xvect[2,]
  x_bis = x*cos(theta) - y*sin(theta)
  y_bis = x*sin(theta) + y*cos(theta)
  (1/10)*(y_bis^2) + x_bis  
}

MetropolisH = function(T_iterations,step,X0,theta,L_k){
  
  chain = matrix(0,nrow=d,ncol=T_iterations)
  
  X_curr = X0
  
  accept = 0
  
  chain[,1] = X_curr
  
  for(i in 2:T_iterations){
    X_cand = (X_curr + step*matrix(rnorm(d,0,1),nrow=d,ncol=1)) / sqrt(1+step^2)
    
    log_r = -(1/2) * (sum(X_cand^2) -  sum(X_curr^2) ) 
    
    u = runif(1)
    
    if(log(u) < log_r){
      if(phi(X_cand,theta) > L_k){
        X_curr = X_cand
        accept = accept + 1
      }
    }
    
    chain[,i] = X_curr
  }
  
  
  acceptation = accept / (T_iterations-1)
  
  return(list(chain=t(chain),taux_accept=acceptation))
  
}

T_iterations  = 20
N = 2000
d = 2

theta_testes = seq(pi/8,2*pi,by=pi/8)
J = 16
tau = 4
#theta = pi/2
essais = 5

N_testes = c(3000)


sapply(N_testes, function(N){
  res_union_plus_duration = sapply(1:essais,function(b){
    start_time <- proc.time()[[3]]
    
    probas = lapply(theta_testes,function(theta){
      
      
      X = matrix(rnorm(d*N,0,1),nrow=d,ncol=N)
      #X0 = X[,1]
      
      p_0 = 0.25
      
      k = 1
      
      phi_X = phi(X,theta)
      
      phi_index = rbind(phi_X,1:N)
      
      phi_index_sort = phi_index[,order(phi_index[1,],decreasing = F)]
      
      L_k = quantile(phi_index_sort[1,],1-p_0)
      
      
      step = 0.3
      
      
      
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
          X_0 = matrix(X_0,nrow=d)
          mec = MetropolisH(T_iterations,step,X_0,theta,L_k)
          mec$chain
        })
        
        #taux_acceptations = res[2,]
        #unlist(taux_acceptations)
        
        #chaines = res[1,]
        
        new_vectors = sapply(1:length_non_accepte,function(i){
          t(chaines[[i]][T_iterations,])
        })
        
        phi_new = phi(new_vectors,theta) 
        
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
          X_0 = matrix(X_0,nrow=d)
          mec = MetropolisH(T_iterations,step,X_0,theta,L_k)
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
    
    
    p_j
    p_j / sum(p_j)
    
    
    Samples_conditionnaly = lapply(1:J,function(i){
      probas[[i]]$X
    })
    
    #View(t(Samples_conditionnaly[[1]])%*%omega )
    
    sum_pj = sum(p_j)
    
    prob_vect = p_j / sum_pj
    
    #r = (N*J)*1/2*T_iterations/100
    
    r = N*J*(1/2)
    
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
    
    s = apply(Xi,2,function(x){
      x = matrix(x,nrow=d)
      calcul = sapply(theta_testes,function(theta){
        phi(x,theta)
      })
      sum(calcul > tau)
    })
    
    res_union = mean( sum_pj / s )
    
    #res_union
    
    end_time <- proc.time()[[3]]
    duration = end_time - start_time
    
    matrix(c(res_union,duration),nrow=2)
    
  })
  
  res_union = res_union_plus_duration[1,]
  duration = res_union_plus_duration[2,]
  
  mean_res_union = mean(res_union)
  sd_res_union = sd(res_union)
  mean_duration = mean(duration)
  
  matrix(c(mean_res_union,sd_res_union,mean_duration),nrow=3)  
  
})




max_phi = function(X)
{
  X = matrix(X,nrow=d)
  max(sapply(theta_testes,phi,Xvect=X))
}


MetropolisH_max = function(T_iterations,step,X0,L_k){
  
  chain = matrix(0,nrow=d,ncol=T_iterations)
  
  X_curr = X0
  
  accept = 0
  
  chain[,1] = X_curr
  
  for(i in 2:T_iterations){
    X_cand = (X_curr + step*matrix(rnorm(d,0,1),nrow=d,ncol=1)) / sqrt(1+step^2)
    
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


#Avec la fonction max directement
essais = 5

N_intuit = c(10000)
#N_intuit = c(100,1000)

sapply(N_intuit,function(N){
  pred_plus_time = sapply(1:essais,function(b){
    start_time <- proc.time()[[3]]  
    
    #vect = omega[,1]  
    
    X = matrix(rnorm(d*N,0,1),nrow=d,ncol=N)
    
    p_0 = 0.25
    
    k = 1
    
    phi_X = apply(X,2,max_phi)
    
    phi_index = rbind(phi_X,1:N)
    
    phi_index_sort = phi_index[,order(phi_index[1,],decreasing = F)]
    
    L_k = quantile(phi_index_sort[1,],1-p_0)
    
    
    step = 0.3
    
    
    
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
        mec = MetropolisH_max(T_iterations,step,X_0,L_k)
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
    
    end_time <- proc.time()[[3]]  
    duration = end_time - start_time
    
    matrix(c(p_hat,duration),nrow=2)
    
  })
  
  pred = pred_plus_time[1,]
  duration = pred_plus_time[2,]
  
  mean_pred = mean(pred)
  
  sd_pred = sd(pred)
  
  mean_dur = mean(duration)
  
  matrix(c(mean_pred,sd_pred,mean_dur),nrow=3)  
})


durations_intuit = c(1.112000e+00,1.343000e+01,5.218600e+01,1.541180e+02)
durations_ALOE = c(3.462000e+00,5.657600e+01,1.165880e+02,2.509520e+02)

sds_intuit = c(4.649702e-07,1.526987e-07,5.788182e-08,6.377317e-08)
sds_ALOE = c(6.350499e-08,1.751369e-08,1.219642e-08 ,6.029530e-09)

df_log_ALOE = data.frame(dur_ALOE=log(durations_ALOE),sds_ALOE=log(sds_ALOE))
df_log_intuit = data.frame(dur_intuit=log(durations_intuit),sds_intuit=log(sds_intuit))

ggplot(df_log_ALOE) + aes(x = dur_ALOE,y=sds_ALOE) + geom_point(color = 'red',size=3) +
  geom_smooth(method = "lm", se = FALSE,col="red") +
  theme_bw() + geom_point(data=df_log_intuit,mapping = aes(x=dur_intuit,y=sds_intuit),color = 'blue',size=3) +
  geom_smooth(data=df_log_intuit,mapping = aes(x=dur_intuit,y=sds_intuit),method = "lm", se = FALSE,col="blue") +
  labs(x = TeX('$\\log_{10}$(computation time)'), 
       y = TeX('$\\log_{10}$(standard deviation)'),
       title =TeX( "Evolution in logarithmic scale of the standard deviation as a function of the computational time"))

