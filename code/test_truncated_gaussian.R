rm(list=objects())
graphics.off()
library(ggplot2)
library(mvtnorm)

library(truncnorm)

x = seq(-5,5,0.1)
density = dtruncnorm(x, a=-Inf, b=Inf, mean = 0, sd = 1)


plot(x,density)




#dim de x
n = 2
r = 1000

x1 = c()
x2 = c()

for (i in 1:r){
  
  z = rnorm(n,0,1)
  
  u = runif(1)
  
  tau = 10
  y = qnorm((u*dnorm(-tau)))
  
  ##on doit avoir t(w)%*%w = 1
  w1 = sqrt(0.25)
  w2 = sqrt(0.75)
  w = matrix(c(w1,w2),ncol=1,nrow = 2)
  x = w%*%y + z - w%*%(t(w)%*%z)
  
  x = -x
  
  x1 = c(x1,x[1,1])
  x2 = c(x2,x[2,1])
}

res = matrix(c(x1,x2),nrow=r,ncol=2,byrow=FALSE)

plot(x1,x2,xlim=c(-20,20),ylim=c(-20,30))
intercept = tau/w2
slope = -w1/w2
abline(a=intercept,b=slope)

qnorm((10**(-12)*pnorm(-10)))

(1e-12*pnorm(-10))


density_conditional = function(x)
{
  denom = pnorm(-tau)
  num = if (t(w)%*%x > tau) 1 else 0
  return((dmvnorm(as.vector(x))*num)/denom)
}

x = matrix(c(5,5),ncol = 1,nrow=2)
density_conditional(x)

1 - pnorm((tau-t(w)%*%x)/(1-t(w)%*%w))



plot(rmvnorm(1000,mean=rep(0,2),sigma=diag(2)),xlim=c(-20,20),ylim=c(-20,30))
points(x1,x2,col='red')



((tau-t(w)%*%x)/(1-t(w)%*%w))




rm(list=objects())
graphics.off()
##First example : circumscribed polygon
library(plotrix)
tau = 1
J = 6
sin_w = sin(2*pi*(1:J)/J)
cos_w = cos(2*pi*(1:J)/J)
omega_transpose = matrix(c(sin_w,cos_w),nrow=J,ncol=2,byrow = F) 
intercept = tau / omega_transpose[,2]
slope = - omega_transpose[,1] / omega_transpose[,2]

par(pty="s") #square plot region

plot(0,0,xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))

mapply(abline, a=intercept,b=slope,col="red")

draw.circle(0,0,tau)

#w1 (en colonne) est bien la direction orthogonale à la droite w1^t * x
lines(c(0,sin(2*pi/J)),c(0,cos(2*pi/J)),col='red')

#comparaison résultats

J = 360
tau = 6

#upper bound de mu
up_bound = exp(-(tau^2)/2)
up_bound

##lower bound 
low_bound = (1- (pi^2)*(tau^2)/(6*J^2)) * up_bound
low_bound

#rapport lower bound / upper bound
rap = (1- (pi^2)*(tau^2)/(6*J^2))
rap

##write ALOE sampling

##Dans notre cas, tous les évènements rares ont autant de chance de se produire

alpha = 1/J ##sert à rien juste pour la compréhension



#dimension
n = 2

sin_w = sin(2*pi*(1:J)/J)
cos_w = cos(2*pi*(1:J)/J)
omega_transpose = matrix(c(sin_w,cos_w),nrow=J,ncol=2,byrow = F) 
taus = 2:8
up_bound = exp(-(taus^2)/2)

r = 100000
set.seed(1)

  
j = sample(1:J,r,replace=T)
  
z = rnorm(n*r,0,1)
  
u = runif(r)
  
y = sapply(taus, function(t){
  qnorm((u*pnorm(-t)))
})

#y = qnorm((u*pnorm(-2)))
  
##on doit avoir t(w)%*%w = 1
w1 = omega_transpose[j,1]
w2 = omega_transpose[j,2]
w = cbind(w1,w2)
wt_z = z[1:r] * w1 + z[(r+1):(r*n)] * w2
z = matrix(z,ncol=2,nrow=r,byrow = F)

#x = (y*w) + z - (wt_z * w)

x = sapply(1:7,function(i){
  (y[,i]*w) + z - (wt_z * w)
})
x = -x
##Les deux colonnes sont mises à la suite

x1 = x[1:r,1] 
x2 = x[(r+1):(2*r),1]
#y



#par(mfrow=c(1,1))
#plot(x1,x2)
#rayon = x1^2 + x2^2
#(rayon > tau^2) ##tous les points simulés sont bien en dehors du cercle de rayon tau
#draw.circle(0,0,tau)

w = omega_transpose

s_x_taus = sapply(1:7, function(k){
  sapply(1:r, function(i){
    wtx = w[,1] * x[i,k] + w[,2] * x[i+r,k]
    sum(wtx > taus[k])
  })
})
  



##Notre estimateur
mu_hat = J*pnorm(-taus[5])*sum(1/s_x_taus[,5]) / r
mu_hat
up_bound[5]


mu_hat/up_bound[5]
#1.016 c'est ok


moy_inv_s_x = sapply(1:7, function(i){
  repartition = matrix(1/(s_x_taus[,i]),nrow=1000,ncol=100)
  apply(repartition,2,mean)
})

mu_hat_vect = sapply(1:7, function(i){
  J*pnorm(-taus[i])*moy_inv_s_x[,i]
})

mu_hat_vect

#mu_hat_vect = J*pnorm(-tau)*moy_inv_s_x

estimate = sapply(1:7,function(i){
  mu_hat_vect[,i] / up_bound[i]
})
estimate
#estimate = mu_hat_vect/up_bound[1]
hist(estimate[,5],breaks=22)

##pvnorm
set.seed(1)
r = 100
mu_hat_mvtnorm = rep(0,r)
tau_vect = rep(tau,J)
sigma = (omega_transpose)%*%t(omega_transpose)
mean = rep(0,J)

tau_vect = sapply(taus, function(x){
  rep(x,J)
})

test = matrix(0,nrow = 360,ncol=r)

mu_mvtnorm = sapply(1:7,function(i){
  apply(test,2,function(x){
    (1-pmvnorm( upper=tau_vect[,i], mean = x, sigma=sigma ))
  })  
})

estimate_mvtnorm = sapply(1:7,function(i){
  mu_mvtnorm[,i] / up_bound[i]
})
 
estimate_mvtnorm

hist(estimate_mvtnorm[,5],breaks=22)

par(mfrow=c(1,2))
hist(estimate[,5],breaks=20)
hist(estimate_mvtnorm[,5],breaks=20)

#Ahn and Kim method

sample = 100

#We sample uniformly at random on the unit sphere
X = rnorm(n*sample)
matrix_x = matrix(X,nrow=sample,ncol=n)
norm_x = apply(matrix_x,1,function(x){
  sqrt(sum(x^2))
})
Theta = (matrix_x / norm_x)

#we compute r(theta)

scalar_product_matrix = omega_transpose%*%t(Theta)

#here tau_j or b_j is constant and mu=0 so the distance min corespond to the greater scalar product
candidate = apply(scalar_product_matrix,2,max)
#ici que des cas positifs car the polytope is bounded in every direction

r_theta = sapply(taus, function(t){
  t/candidate
})

r_theta_squared = r_theta^2

proba_R_Rtheta = 1 - pchisq(r_theta_squared,n)

mean_pr_R_Rtheta = apply(proba_R_Rtheta,2,mean)

#rapport on the upper bound (considered as exact)
mean_pr_R_Rtheta / up_bound

estimate_AK = proba_R_Rtheta%*%diag(1/up_bound)

sd_error = apply(proba_R_Rtheta,2,sd)
sd_error

##asymptotically efficient/optimal
estimateur_proba_carre = proba_R_Rtheta^2

mean_pr_carre_R_Rtheta = apply(estimateur_proba_carre,2,mean)

num = log(mean_pr_carre_R_Rtheta)

den = log(mean_pr_R_Rtheta^2)

ratio = num/den
ratio

##mean square relative error

msre_AK = apply(estimate_AK,2,function(x){
  mean((x - 1)^2)
})
msre_AK

##ALOE
msre_aloe = apply(estimate,2,function(x){
  mean((x - 1)^2)
})
#msre_aloe = mean( (estimate - 1)^2 ) 
msre_aloe
##pmvnorm
msre_mvtnorm = apply(estimate_mvtnorm,2,function(x){
  mean((x - 1)^2)
})
#msre_mvtnorm = mean( (estimate_mvtnorm - 1)^2 ) 
msre_mvtnorm

rapport = msre_mvtnorm / msre_aloe
rapport

##upper bound on rmse, en prenant mu = upper bound de mu = exp(-tau^2 / 2)

##on avait pris n égal 1000 pour construire un seul estimateur
n = 1000
#valeur donnée pour tau=6
msre_bound = ((J*pnorm(-6)/up_bound[5]) - 1) / n
msre_bound
# 0.02232055
##This is in line with the figure given in the paper. So to compute only one estimator, he indeed used 
#n = 1000.
msre_bound/msre_aloe[5] 
##pas vraiment 20..

##On change un peu l'exemple pour tuer la symétrie

PrimeNumber <- function(n){
  #Eratosthenes 
  #Return all prime numbers up to n (based on the sieve of Eratosthenes)
  if (n >= 2) {
    sieve <- seq(2, n)
    primes <- c()
    
    for (i in seq(2, n)) {
      if (any(sieve == i)) {
        primes <- c(primes, i)
        sieve <- c(sieve[(sieve %% i) != 0], i)
      }
    }
    return(primes)
  } else {
    stop("Input value of n should be at least 2.")
  }
}

#we generate all the prime numbers among integers up to 360
PrimeNumber(360)

angles = 2*pi*PrimeNumber(360) / 360
J = length(angles)
dim = 2
#J = 72

sin_w = sin(angles)
cos_w = cos(angles)
omega_transpose = matrix(c(sin_w,cos_w),nrow=J,ncol=2,byrow = F) 

r = 100000
set.seed(1)


j = sample(1:J,r,replace=T)

z = rnorm(dim*r,0,1)

u = runif(r)

y = sapply(taus, function(t){
  qnorm((u*pnorm(-t)))
})

#y = qnorm((u*pnorm(-2)))

##on doit avoir t(w)%*%w = 1
w1 = omega_transpose[j,1]
w2 = omega_transpose[j,2]
w = cbind(w1,w2)
wt_z = z[1:r] * w1 + z[(r+1):(r*dim)] * w2
z = matrix(z,ncol=2,nrow=r,byrow = F)

#x = (y*w) + z - (wt_z * w)

x = sapply(1:7,function(i){
  (y[,i]*w) + z - (wt_z * w)
})
x = -x
##Les deux colonnes sont mises à la suite

x1 = x[1:r,1] 
x2 = x[(r+1):(2*r),1]
#y



par(mfrow=c(1,1))
plot(x1,x2)
rayon = x1^2 + x2^2
(rayon > tau^2) ##tous les points simulés sont bien en dehors du cercle de rayon tau
draw.circle(0,0,tau)

w = omega_transpose

s_x_taus = sapply(1:7, function(k){
  sapply(1:r, function(i){
    wtx = w[,1] * x[i,k] + w[,2] * x[i+r,k]
    sum(wtx > taus[k])
  })
})




moy_inv_s_x = sapply(1:7, function(i){
  repartition = matrix(1/(s_x_taus[,i]),nrow=1000,ncol=100)
  apply(repartition,2,mean)
})

mu_hat_vect = sapply(1:7, function(i){
  J*pnorm(-taus[i])*moy_inv_s_x[,i]
})

mu_hat_vect

#mu_hat_vect = J*pnorm(-tau)*moy_inv_s_x

estimate = sapply(1:7,function(i){
  mu_hat_vect[,i] / up_bound[i]
})
estimate

##ALOE
msre_aloe = apply(estimate,2,function(x){
  mean((x - 1)^2)
})
#msre_aloe = mean( (estimate - 1)^2 ) 
msre_aloe

n = 1000
#valeur donnée pour tau=6
msre_bound = ((J*pnorm(-6)/up_bound[5]) - 1) / n
msre_bound

##variance de mu_hat/mu (avec mu choisi comme up_bound) 
var(estimate[,5])



mean(mu_hat_vect[,5]) / up_bound[5] 
##la moyenne des estimateurs donne comme estimation 0.97 de l'upper bound

##Ahn and Kim estimator

samples = 100

#We sample uniformly at random on the unit sphere
X = rnorm(dim*samples)
matrix_x = matrix(X,nrow=samples,ncol=dim)
norm_x = apply(matrix_x,1,function(x){
  sqrt(sum(x^2))
})
Theta = (matrix_x / norm_x)

#we compute r(theta)

scalar_product_matrix = omega_transpose%*%t(Theta)

#here tau_j or b_j is constant and mu=0 so the distance min corespond to the greater scalar product
candidate = apply(scalar_product_matrix,2,max)
#ici que des cas positifs car the polytope is bounded in every direction

r_theta = sapply(taus, function(t){
  t/candidate
})

r_theta_squared = r_theta^2

proba_R_Rtheta = 1 - pchisq(r_theta_squared,dim)

mean_pr_R_Rtheta = apply(proba_R_Rtheta,2,mean)
mean_pr_R_Rtheta

#rapport on the upper bound (considered as exact)
mean_pr_R_Rtheta / up_bound

estimate_AK = proba_R_Rtheta%*%diag(1/up_bound)

sd_error = apply(proba_R_Rtheta,2,sd)
sd_error

##asymptotically efficient/optimal
estimateur_proba_carre = proba_R_Rtheta^2

mean_pr_carre_R_Rtheta = apply(estimateur_proba_carre,2,mean)

num = log(mean_pr_carre_R_Rtheta)

den = log(mean_pr_R_Rtheta^2)

ratio = num/den
ratio

##pvnorm
set.seed(1)
r = 100
mu_hat_mvtnorm = rep(0,r)
sigma = (omega_transpose)%*%t(omega_transpose)
mean = rep(0,J)

tau_vect = sapply(taus, function(x){
  rep(x,J)
})

test = matrix(0,nrow = J,ncol=r)

mu_mvtnorm = sapply(1:7,function(i){
  apply(test,2,function(x){
    (1-pmvnorm( upper=tau_vect[,i], mean = x, sigma=sigma ))
  })  
})

estimate_mvtnorm = sapply(1:7,function(i){
  mu_mvtnorm[,i] / up_bound[i]
})

estimate_mvtnorm

msre_mvtnorm = apply(estimate_mvtnorm,2,function(x){
  mean((x - 1)^2)
})
#msre_mvtnorm = mean( (estimate_mvtnorm - 1)^2 ) 
msre_mvtnorm

##variance de mu_hat/mu (avec mu choisi comme up_bound) dans le cas mvtnorm 
var(estimate_mvtnorm[,5])


####Example High dimensional half-spaces

##on construit 200 problèmes et pour chaque on calcul un estimateur de mu


m = 200

d_list = c(20,50,100,200,500)

d = sample(d_list,m,replace = T)

J_list = c(1/2,1,2)

J = sample(J_list,m,replace=T)

J = J*d

#set.seed(1)

##generate random unit vectors
##d n'est jamais impair

w = sapply(1:m,function(i){
  x = rnorm(d[i]*J[i])
  matrix_x = matrix(x,nrow=d[i],ncol=J[i])
  norme = apply(matrix_x,2,function(x){
    sqrt(sum(x^2))
  })
  ##to multiply each column by a value, we use diag matix with the vector and matrix multiplication
  wj = matrix_x%*%diag(1/norme)
})

log_10_union_bound = -runif(m,4,8)

union_bound = 10^(log_10_union_bound)

tau = -qnorm((union_bound)/ J)


##Ahn and Kim method

samples = 1000

Theta_list = sapply(d,function(di){
  #We sample uniformly at random on the unit sphere
  X = rnorm(di*samples)
  matrix_x = matrix(X,nrow=samples,ncol=di)
  norm_x = apply(matrix_x,1,function(x){
    sqrt(sum(x^2))
  })
  Theta = (matrix_x / norm_x)
})


#we compute r(theta)
r_theta_list = sapply(1:m,function(i){
  scalar_product_matrix = Theta_list[[i]]%*%w[[i]]
})

#here tau_j or b_j is constant and mu=0 so the distance min corespond to the greater scalar product
candidate_matrix = sapply(1:m,function(i){
  candidate = apply(r_theta_list[[i]],1,max)
})

sum(candidate_matrix < 0)
#4 , in the 200 problems with each 100 directions sampled, there are 4 directions which are not bounded by any 
#half space in the problem

which(candidate_matrix<0)

candidate_matrix[11416]

T = sapply(1:m,function(i){
  r_th = tau[i]/candidate_matrix[,i]
  r_th_squared = ifelse(r_th  > 0,r_th^2,0)
  prob_R_rtheta = ifelse(r_th_squared  > 0,pchisq(r_th_squared,d[i],lower.tail = F),0)  
})


mean_pr_R_rtheta = apply(T,2,mean)

lower_bound = pnorm(-tau)
upper_bound = J*pnorm(-tau)

estimate_AK = T%*%diag(1/lower_bound)
results = apply(estimate_AK,2,mean)
results

results_up_bound = apply(T%*%diag(1/upper_bound),2,mean)

results_up_bound[which(results_up_bound > 2)]

hist(log10(results),breaks=20)



#ALOE
n = 1000 
total_sample = m*n


##on est encore dans le cas où les évènements rares sont tirés de manière uniforme
j = sapply(J,function(x){
  sample(1:x,n,replace=T)
})
#pour chaque problème une colonne de longueur 1000 de j

z = sapply(d,function(x){
  rnorm(n*x,0,1)
})
#z, liste de de longueur 200 (1 élément par problème) , n=1000 vecteurs de dimension d

y = sapply(tau, function(t){
  u = runif(n)
  qnorm((u*pnorm(-t)))
})
##matrice de y, une colonne pour chaque problème

x = sapply(1:m,function(i){
  omega = w[[i]][,j[,i]]  
  matrix_z = matrix(z[[i]],nrow=n,ncol=d[i])
  wtz_vect = diag(matrix_z%*%omega)
  -(omega%*%diag(y[,i]) + t(matrix_z) - omega%*%diag(wtz_vect))
})
##liste de matrices, pour chaque matrice la colonne est l'un vecteur de dimension d

##Calcul de S(x)



s_x = sapply(1:m, function(i){
  matrix_wtx = t(x[[i]])%*%w[[i]]
  s = apply(matrix_wtx,1,function(x){
    sum(x > tau[i])
  })
})

mean(apply(s_x,2,sum) == 1000) * 100
##81.5 % of the simulations (of the problems) had no intersections among 1000 trials. 


mean_inv_s_x = apply(s_x,2,function(x){
  mean(1/x)
})

#union_bound = J*pnorm(-tau)

mu_hat = union_bound*mean_inv_s_x

rate_mu_hat_union_bound =  mean_inv_s_x

plot(log10(union_bound),rate_mu_hat_union_bound)
#axis in x is logarithmic

##pvmnorm

set.seed(1)
sigma = sapply(1:m,function(i){
  t(w[[i]])%*%(w[[i]])
})

mean = sapply(J,function(x){
  rep(0,x)
})

tau_vect = sapply(1:m,function(i){
  rep(tau[i],J[i])
})

#pmvnorm( upper=tau_vect, mean = mean, sigma=sigma )

mu_hat_mvtnorm = sapply(1:m,function(i){
  (1-pmvnorm( upper=tau_vect[[i]], mean = mean[[i]], sigma=sigma[[i]] ))
})

estimate_rate_mvtnorm = mu_hat_mvtnorm / union_bound

library(latex2exp) # permet d'ajouter des equations latex aux figures

plot(log10(union_bound),estimate_rate_mvtnorm,col='red',pch=19,
main="Results of 200 estimates of the μ for varying high dimensional problems with nearly independent events",
xlab=TeX("$log_{10}(\\bar{\\mu})$"),ylab="Estimate/Bound")
points(log10(union_bound),rate_mu_hat_union_bound,col='blue')

legend(
  "topright", 
  legend=c("ALOE", "pmvnorm"), 
  bty='o', 
  pch=c(1,19),
  lty=c(0,0),
  col=c('blue','red')
)


hist(log10(rate_mu_hat_union_bound),breaks=20)



###Specific example to illustrate asymptotically efficient estimator

library(latex2exp) # permet d'ajouter des equations latex aux figures
library(ggplot2)

d = 20
J = 2000

x = rnorm(d*J)
matrix_x = matrix(x,nrow=d,ncol=J)
norme = apply(matrix_x,2,function(x){
  sqrt(sum(x^2))
})
##to multiply each column by a value, we use diag matix with the vector and matrix multiplication
omega = matrix_x%*%diag(1/norme)

#log_10_union_bound = -seq(160,1,by=-3)
log_10_union_bound = -seq(60,1,by=-1)

union_bound = 10^(log_10_union_bound)

tau = -qnorm((union_bound)/ J)
tau
l = length(tau)
#tau = 1:20


samples = 1000


#We sample uniformly at random on the unit sphere
X = rnorm(d*samples)
matrix_x = matrix(X,nrow=samples,ncol=d)
norm_x = apply(matrix_x,1,function(x){
  sqrt(sum(x^2))
})
Theta = (matrix_x / norm_x)

#we compute r(theta)

scalar_product_matrix = Theta%*%omega

#here tau_j or b_j is constant and mu=0 so the distance min corespond to the greater scalar product
candidate = apply(scalar_product_matrix,1,max)

sum(candidate < 0)
  
estimateur_one_sample = sapply(tau,function(t){
  r_th = t/candidate
  r_th_squared = ifelse(r_th  > 0,r_th^2,0)
  prob_R_rtheta = ifelse(r_th_squared  > 0,pchisq(r_th_squared,d,lower.tail = F),0) 
})

  
num = log(apply(estimateur_one_sample^2,2,mean))

den = log( (apply(estimateur_one_sample,2,mean))^2 )

ratio = exp(log(-num)-log(-den))
ratio

df_u_bound = data.frame(log_10=log_10_union_bound,ratio=ratio)
df_tau = data.frame(tau=tau,ratio=ratio)

plt_ubound = ggplot(df_u_bound) +
      aes(x = log_10_union_bound) +
      geom_point(mapping = aes(y=ratio)) +
      geom_hline(yintercept=1, colour='blue') +
      theme_bw() +
      labs(x = TeX('$\\log_{10}(\\bar{\\mu})$'), 
           y = TeX('$l(\\tau)$'),
           title = TeX("Evolution of the ratio $l(\\tau)$ as a function of $\\log_{10}(\\bar{\\mu})$"))

plt_ubound

plt_tau = ggplot(df_tau) +
  aes(x = tau) +
  geom_point(mapping = aes(y=ratio)) +
  geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\tau$'), 
       y = TeX('$l(\\tau)$'),
       title = TeX("Evolution of the ratio $l(\\tau)$ as a function of $\\tau$"))

plt_tau


plot(log_10_union_bound,ratio,xlab=TeX('$\\log_{10}(\\bar{\\mu})$'),ylab = TeX('$\\frac{\\log(\mathbb{E}[L(\\tau)^2])}{\\log(\mathbb{E}[L(\\tau)]^2)}$'))

plot(tau,ratio)


##strongly efficient estimator for the ALOE method

##on est encore dans le cas où les évènements rares sont tirés de manière uniforme

j = sample(1:J,samples,replace=TRUE)

z = rnorm(samples*d,0,1)

y = sapply(tau, function(t){
  u = runif(samples)
  qnorm((u*pnorm(-t)))
})

omega_rand = omega[,j]  
matrix_z = matrix(z,nrow=samples,ncol=d)
wtz_vect = diag(matrix_z%*%omega_rand)

x = lapply(1:l,function(i){
  -(omega_rand%*%diag(y[,i]) + t(matrix_z) - omega_rand%*%diag(wtz_vect))
})

##Calcul de S(x)
s_x = sapply(1:l, function(i){
  matrix_wtx = t(x[[i]])%*%omega
  s = apply(matrix_wtx,1,function(x){
    sum(x > tau[i])
  })
})


expectation = sapply(1:l, function(i){
  mean( union_bound[i] / s_x[,i] )
})

expectation_square = sapply(1:l, function(i){
  mean( union_bound[i]^2 / s_x[,i]^2 )
})

expect_inv_s_x_square = apply(s_x,2,function(x){
  mean( (1/x)^2 )
})

#union_bound = J*pnorm(-tau)


expect_inv_s_x = apply(s_x,2,function(x){
  mean( 1/x )
})

ratio = (expect_inv_s_x_square - expect_inv_s_x^2) / (expect_inv_s_x^2)
ratio

df_aloe_ratio = data.frame(tau=tau,ratio=ratio)
df_aloe_ratio_log10 = data.frame(log10 = log_10_union_bound,ratio=ratio)

plt_aloe = ggplot(df_aloe_ratio) +
  aes(x = tau) +
  geom_point(mapping = aes(y=ratio)) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\tau$'), 
       y = TeX('$r(\\tau)$'),
       title = TeX("Evolution of the ratio $r(\\tau)$ as a function of $\\tau$"))

plt_aloe


plt_log10 = ggplot(df_aloe_ratio_log10) +
  aes(x = log10) +
  geom_point(mapping = aes(y=ratio)) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\log_{10}(\\bar{\\mu})$'), 
       y = TeX('$r(\\tau)$'),
       title = TeX("Evolution of the ratio $r(\\tau)$ as a function of $\\log_{10}(\\bar{\\mu})$"))

plt_log10



##ratio de 1 car en grande dimension et peu de J par rapport à cette grande dimension (donc évènements rares ne se
##superposent quasiment pas). (mu presque égal à l'union bound)


#théoriquement comme construit avec n=1 à chaque fois (pour les 1000 samples) la borne théorique est J+n-1 / n
(J-1)

###rapport des moments pour la méthode de Ahn and Kim

second_moment_a_k = apply(estimateur_one_sample^2,2,mean)

first_moment_a_k = (apply(estimateur_one_sample,2,mean))

strong_efficient_directional = (second_moment_a_k - (first_moment_a_k)^2) / (first_moment_a_k)^2
strong_efficient_directional

df_dir_s_ef_tau = data.frame(tau=tau,res=strong_efficient_directional)
df_dir_s_ef_mu = data.frame(log10=log_10_union_bound,res=strong_efficient_directional)


plt_ds_tau = ggplot(df_dir_s_ef_tau) +
  aes(x = tau) +
  geom_point(mapping = aes(y=res)) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\tau$'), 
       y = TeX('$r_d(\\tau)$'),
       title = TeX("Evolution of the ratio $r_d(\\tau)$ as a function of $\\tau$"))

plt_ds_tau


plt_ds_log10 = ggplot(df_dir_s_ef_mu) +
  aes(x = log10) +
  geom_point(mapping = aes(y=res)) +
  #geom_hline(yintercept=1, colour='blue') +
  theme_bw() +
  labs(x = TeX('$\\log_{10}(\\bar{\\mu})$'), 
       y = TeX('$r_d(\\tau)$'),
       title = TeX("Evolution of the ratio $r_d(\\tau)$ as a function of $\\log_{10}(\\bar{\\mu})$"))

plt_ds_log10

