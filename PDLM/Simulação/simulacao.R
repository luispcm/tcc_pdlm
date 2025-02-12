library(tidyverse)
library(MASS)
library(Matrix)
library(corpcor)
library(mvtnorm)
library(scales)
library(ggthemes)

################ Gerando os dados (Curva 1) #############
set.seed(80)

x0 = 60
x=c()

#desvio padrão do erro das regressoras:
tau2_raiz = sqrt(10)

#passo 1
x[1] = x0 + rnorm(1,0,tau2_raiz)

n = 365 #quantidade de dados
for (i in 2:n){
  
  x[i] = x[i-1] + rnorm(1,0,tau2_raiz)
  
}

plot(x)
#definindo os betas e o intercepto:

#vetor beta com dimensão (Lags)

lags = 16
q = lags + 1 #inclui o beta0

#montando a matriz de regressoras defasadas

X = matrix(nrow = n,
           ncol = q)

X[,1] = x #lag = 0

for (i in 2:q){
  
  X[,i] = dplyr::lag(x, i-1)
  
}

#X efetivo (sem os NA's)::

X = X[-c(1:(q-1)),]

######################## Curva 1 (d = 3) ####################

#vetor de etas verdadeiros:

eta_ver = c(0.95, 0.59, 0.1, -0.01)/100

M_aux = matrix(nrow = q-2,
               ncol = 4)

for (i in 2:(q-1)){
  
  M_aux[i-1,] =  c(1, i, i^2, i^3)
  
}

beta_ver = c(eta_ver[1],
             sum(eta_ver))

for (i in 1:(q-2)){
  beta_ver = c(beta_ver,eta_ver%*%M_aux[i,])
}

plot(x = 1:q, y = beta_ver)
abline(h= 0)

#contruindo o St,i

#grau do polinômio d = 3

v = c(1:(q-1)) #vetor que multiplica cada matriz

St0 = apply(X,1,sum)
St1 = X[,-1]%*%v^1
St2 = X[,-1]%*%v^2
St3 = X[,-1]%*%v^3

############## Modelo Linear Dinâmico com Defasagem Polinomial #########

time = nrow(X) #tamanho da série temporal
n = 5       # Dimensão do vetor dos parâmetros a cada tempo

#vetor F

Ft = matrix(rep(1),n,time)

Ft[2,] = St0
Ft[3,] = St1
Ft[4,] = St2
Ft[5,] = St3

#Bloco da G para o nível

G1 = diag(rep(1,1))

#Bloco da G para as regressoras

G2 = diag(rep(1,n-1)) #n-1 pois tiramos o nível

#matriz G

Gt = bdiag(G1, G2)
Gt = array(as.matrix(Gt), dim = c(n,n,time))
Gt

#### gerando o Y

alpha_1 = 20

V_raiz = sqrt(5)

W = diag(c(1, rep(0, n-1))) #matrix de variâncial de evolução

theta_t = matrix(c(alpha_1,eta_ver), nrow = n, ncol = time)

Y = t(Ft[,1])%*%theta_t[,1] + rnorm(1, 0, V_raiz)

for (i in 2:time){
  
  #equação de estados
  theta_t[,i] = Gt[,,i]%*%theta_t[,i-1] + mvrnorm(1, rep(0, n), W)
  
  #equação de evolução
  
  Y[i] = t(Ft[,i])%*%theta_t[,i] + rnorm(1, 0, V_raiz)
  
}

alpha = theta_t[1,]
plot.ts(theta_t[1,])
abline(h=0)

plot.ts(Y)

##### Estimando a série

# Padronizando as variáveis

Ft_padronizado = apply(Ft,1,sd)
Ft_padronizado[1] = 1
 # a padronização ajuda a evitar problemas numéricos na estimação
Ft = Ft/Ft_padronizado

#Criando os componentes

at = array(NA, dim=c(n,1,time))
Rt = array(NA, dim=c(n,n,time))

mt = array(NA, dim=c(n,1,time))
Ct = array(NA, dim=c(n,n,time))

ft = rep(NA, time)
Qt = rep(NA, time)

et = rep(NA, time)
At = array(NA, dim=c(n,1,time))

St = rep(NA,time)
nt = rep(NA,time)

##### Fator de Desconto

#Bloco do nível

d1 = 0.99
fd1 = 1/d1

#Bloco da regressora

d2 = 1
fd2 = 1/d2

#Juntando tudo

D_inv = diag(c(fd1, rep(fd2, n-1)))
D_inv

####### Filtro de Kalman

#Passo t  = 0

m0 = c(0, rep(0, n-1))
C0 = diag(c(100, rep(100,n-1)), n, n)

n0 = 3
S0 = 0.5

curve(dgamma(x,n0/2,n0*S0/2), from =0 , to = 20)
abline(v = 1/V_raiz, col = "red")

### Passo t=1

# Priori em t=1
at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%t(Gt[,,1])) + (D_inv-diag(1, n))*(Gt[,,1]%*%C0%*%t(Gt[,,1]))

# Previsão 1 passo-a-frente
ft[1] = t(Ft[,1])%*%at[,,1]
Qt[1] = (t(Ft[,1])%*%Rt[,,1]%*%Ft[,1]) + S0

# Posteriori em t=1
At[,,1] = (Rt[,,1]%*%Ft[,1])*(1/Qt[1])
et[1] = Y[1] - ft[1]

#Equações de phi
nt[1] = n0 + 1
St[1] = S0 + (S0/nt[1])*(((et[1])^2/Qt[1])-1)

mt[,,1] = at[,,1] + (At[,,1]*et[1])
Ct[,,1] = (St[1]/S0)*(Rt[,,1] - At[,,1]%*%t(At[,,1])*Qt[1])


for(t in 2:time){ # Passo 2 até o passo TIME
  
  # Priori em t
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] = Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]) + (D_inv-diag(1, n))*(Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]))
  
  # Previsão 1 passo-a-frente
  ft[t] = t(Ft[,t])%*%at[,,t]
  Qt[t] = (t(Ft[,t])%*%Rt[,,t]%*%Ft[,t]) + St[t-1]
  
  # Posteriori em t
  At[,,t] = (Rt[,,t]%*%Ft[,t])%*%(1/Qt[t])
  et[t] = Y[t] - ft[t]
  
  #Equações de phi
  nt[t] = nt[t-1] + 1
  St[t] = St[t-1] + ((St[t-1]/nt[t])*(((et[t])^2/Qt[t])-1))
  
  mt[,,t] = at[,,t] + (At[,,t]*et[t])
  Ct[,,t] = St[t]/St[t-1]*(Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[t])
  
}

#Analisando

plot.ts(Y)
lines(ft, col = "red")
plot.ts(sqrt(Qt[10:t]))

plot.ts(alpha)
lines(mt[1,,], col = "red")

plot(mt[c(2:n),,time]/Ft_padronizado[-1], col = "red")
points(eta_ver, col = "black")


##### Suavização

#definindo Bt:
Bt = array(NA, dim=c(n,n,time))

mts = mt 
Cts = Ct

for (k in 1:(time-1)){
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time-k+1])
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mts[,,time-k+1] - at[,,time-k+1])
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Cts[,,time-k+1] - Rt[,,time-k+1])%*%t(Bt[,,time-k])
  
  if(is.positive.definite(Cts[,,time-k]) == F){Cts[,,time-k] = make.positive.definite(Cts[,,time-k])}
}

plot.ts(alpha)
lines(mt[1,,1:time], col = "blue")
lines(mt[1,, 1:time] + 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mt[1,, 1:time] - 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mts[1,, 1:time], col = "red")
lines(mts[1,, 1:time] + 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)
lines(mts[1,, 1:time] - 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)

#Amotrando thetha

amostra_theta = array(NA, dim = c(1000, n, time))

#no tempo t = 1
amostra_theta[,,1] = rmvt(1000, delta = mts[,,1], sigma = St[1]/S0*Cts[,,1], df = nt[1])

#no tempo t = 2 em diante
for(i in 2:time){
  if(isSymmetric(Cts[,,i]) == F){
    upper_tri <- upper.tri(Cts[,,i])
    # Tornando a matriz simétrica
    Cts[,,i] <- Cts[,,i] * upper_tri + t(Cts[,,i]) * (1 - upper_tri)
  }
  
  amostra_theta[,,i] = rmvt(1000, delta = mts[,,i], sigma = St[i]/St[i-1]*Cts[,,i], df = nt[i])
}

# retomando as escalas

for(i in 1:time){
  amostra_theta[,,i] = t(t(amostra_theta[,,i])/Ft_padronizado) # retirando o efeito compensatório feito pelo modelo devido à padronização feita anteriormente
}

#amostra do vetor eta

amostra_eta = amostra_theta[,c(2:n),time]

plot(apply(amostra_eta, FUN = quantile, MARGIN = 2, probs = 0.5))
points(eta_ver, col = "red")

#amostra do vetor beta (nas análises seguintes eu pego sempre o último vetor estimado, pois como ele é constante, eu devo pegar o mais atualizado)

#matriz auxiliar
m = matrix(1, nrow = q,
           ncol = n-1)

m[1,] = c(1,0,0,0)

for (i in 2:(n-1)){
  for (j in 1:(q-1)){
    
    m[j+1,i] = j^(i-1)
    
  }
}

amostra_beta = amostra_eta%*%t(m)

#y estimado
mu = phi = c()

for (i in 1:time){
  mu[i] = t(Ft[,i])%*%mts[,,i]
  phi[i] = qgamma(0.5, nt[i]/2, nt[i]*St[i]/2)
}

y_estimado = qnorm(0.5,mu,sqrt(1/phi))
ic_inf_y = qnorm(0.025,mu,sqrt(1/phi))
ic_sup_y = qnorm(0.975,mu,sqrt(1/phi))


############### Intervalo de Credibilidade de 95% ##########

#Pego o beta mais atualizado, pois ele é constante ao longo do tempo

ic_inf_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.025)
beta_estimado = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.5)
ic_sup_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.975)

#### Y estimado

ggplot(NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y, color = "Dados"), size = 2) +
  geom_line(aes(y = y_estimado, color = "Estimativas"), size = 2) +
  geom_ribbon(aes(ymin = ic_inf_y, ymax = ic_sup_y,
                  fill = "IC 95%"),
              alpha = 0.3) +
  geom_hline(yintercept = 0) +
  labs(x = "Tempo",
       y = "y",
       colour = "") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 45),
        axis.text = element_text(size = 45,
                                 angle = 0),
        legend.text=element_text(size= 30),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_text(size = 45)) +
  scale_color_manual(values = c("Dados" = "black", "Estimativas" = "red"),
                     breaks = c("Estimativas", "Dados")) +
  scale_fill_manual(values = c("IC 95%" = "blue")) +
  scale_y_continuous(limits = c(0.9*min(Y), 1.1*max(Y)))

############ Gráfico dos betas estimados #########
ggplot(NULL,
       aes(x = 0:(q-1))) +
  geom_point(aes(y = beta_estimado), size = 8) +
  geom_errorbar(aes(ymin = ic_inf_beta, ymax = ic_sup_beta),
                width = 0.5, size = 2) +
  geom_point(aes(y = beta_ver, colour = "red"), size = 8, show.legend = F) +
  labs(x = "Lags",
       y = bquote(beta),
       colour = "") +
  theme_hc(base_size = 1) +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        legend.position = NULL,
        panel.grid = element_line(linewidth = 1)) +
  scale_x_continuous(breaks = seq(0, (q-1), by = 1)) +
  scale_y_continuous(breaks = pretty(range(ic_inf_beta, ic_sup_beta)),
                     labels = label_number(decimal.mark = ","))



################ Gerando os dados (Curva 2, d=2) #############

#vetor de etas verdadeiros:

d = 2

eta_ver = c(2, -0.3, 0.01)/10

M_aux = matrix(nrow = q-2,
               ncol = d+1)

for (i in 2:(q-1)){
  
  M_aux[i-1,] =  c(1, i, i^2)
  
}

beta_ver = c(eta_ver[1],
             sum(eta_ver))

for (i in 1:(q-2)){
  beta_ver = c(beta_ver,eta_ver%*%M_aux[i,])
}

plot(x = 1:q, y = beta_ver)
abline(h= 0)

#contruindo o St,i

#grau do polinômio d = 2

v = c(1:(q-1)) #vetor que multiplica cada matriz

St0 = apply(X,1,sum)
St1 = X[,-1]%*%v^1
St2 = X[,-1]%*%v^2
#St3 = X[,-1]%*%v^3

############## Modelo Linear Dinâmico com Defasagem Polinomial #########

time = nrow(X) #tamanho da série temporal
n = 4        # Dimensão do vetor dos parâmetros a cada tempo

#vetor F

Ft = matrix(rep(1),n,time)

Ft[2,] = St0
Ft[3,] = St1
Ft[4,] = St2
#Ft[5,] = St3

#Bloco da G para o nível

G1 = diag(rep(1,1))

#Bloco da G para as regressoras

G2 = diag(rep(1,n-1)) #n-1 pois tiramos o nível

#matriz G

Gt = bdiag(G1, G2)
Gt = array(as.matrix(Gt), dim = c(n,n,time))
Gt

#### gerando o Y

V_raiz = sqrt(5)

W = diag(c(1, rep(0, n-1))) #matrix de variâncial de evolução

alpha_1 = 20

theta_t = matrix(c(alpha_1,eta_ver), nrow = n, ncol = time)

Y = t(Ft[,1])%*%theta_t[,1] + rnorm(1, 0, V_raiz)

for (i in 2:time){
  
  #equação de estados
  theta_t[,i] = Gt[,,i]%*%theta_t[,i-1] + mvrnorm(1, rep(0, n), W)
  
  #equação de evolução
  
  Y[i] = t(Ft[,i])%*%theta_t[,i] + rnorm(1, 0, V_raiz)
  
}

alpha = theta_t[1,]
plot.ts(theta_t[1,])
abline(h=0)

plot.ts(Y)


##### Estimando a série

# Padronizando as variáveis

Ft_padronizado = apply(Ft,1,sd)
Ft_padronizado[1] = 1

Ft = Ft/Ft_padronizado

#Criando os componentes

at = array(NA, dim=c(n,1,time))
Rt = array(NA, dim=c(n,n,time))

mt = array(NA, dim=c(n,1,time))
Ct = array(NA, dim=c(n,n,time))

ft = rep(NA, time)
Qt = rep(NA, time)

et = rep(NA, time)
At = array(NA, dim=c(n,1,time))

St = rep(NA,time)
nt = rep(NA,time)

##### Fator de Desconto

#Bloco do nível

d1 = 0.99
fd1 = 1/(d1)

#Bloco da regressora

d2 = 1
fd2 = 1/(d2)

#Juntando tudo

D_inv = diag(c(fd1, rep(fd2, n-1)))
D_inv

####### Filtro de Kalman

#Passo t  = 0

m0 = c(0, rep(0, n-1))
C0 = diag(c(100, rep(100,n-1)), n, n)

n0 = 3
S0 = 0.5

curve(dgamma(x,n0/2,n0*S0/2), from =0 , to = 20)
abline(v = 1/V_raiz, col = "red")

#Passo t = 1

at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%t(Gt[,,1])) + (D_inv-diag(1, n))*(Gt[,,1]%*%C0%*%t(Gt[,,1]))

ft[1] = t(Ft[,1])%*%at[,,1]
Qt[1] = t(Ft[,1])%*%Rt[,,1]%*%Ft[,1] + S0

et[1] = Y[1] - ft[1]
At[,,1] = Rt[,,1]%*%Ft[,1]/Qt[1]

nt[1] = n0 + 1
St[1] = S0 + S0/nt[1]*(et[1]^2/Qt[1] - 1)

mt[,,1] = at[,,1] + At[,,1]*et[1]
Ct[,,1] = St[1]/S0*(Rt[,,1] - At[,,1]%*%t(At[,,1])*Qt[1])

#Passo t de 2 a time

for (t in 2:time){
  
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] = Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]) + (D_inv-diag(1, n))*(Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]))
  
  ft[t] = t(Ft[,t])%*%at[,,t]
  Qt[t] = t(Ft[,t])%*%Rt[,,t]%*%Ft[,t] + St[t-1]
  
  et[t] = Y[t] - ft[t]
  At[,,t] = Rt[,,t]%*%Ft[,t]/Qt[t]
  
  nt[t] = nt[t-1] + 1
  St[t] = St[t-1] + St[t-1]/nt[t]*(et[t]^2/Qt[t] - 1)
  
  mt[,,t] = at[,,t] + At[,,t]*et[t]
  Ct[,,t] = St[t]/St[t-1]*(Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[t])
  
}

#Analisando

plot.ts(Y)
lines(ft, col = "red")

plot.ts(alpha)
lines(mt[1,,], col = "red")

plot(mt[c(2:n),,time]/Ft_padronizado[-1], col = "red")
points(eta_ver, col = "black")


##### Suavização

#definindo Bt:
Bt = array(NA, dim=c(n,n,time))

mts = mt 
Cts = Ct

for (k in 1:(time-1)){
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time-k+1])
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mts[,,time-k+1] - at[,,time-k+1])
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Cts[,,time-k+1] - Rt[,,time-k+1])%*%t(Bt[,,time-k])
  
  if(is.positive.definite(Cts[,,time-k]) == F){Cts[,,time-k] = make.positive.definite(Cts[,,time-k])}
}

plot.ts(alpha)
lines(mt[1,,1:time], col = "blue")
lines(mt[1,, 1:time] + 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mt[1,, 1:time] - 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mts[1,, 1:time], col = "red")
lines(mts[1,, 1:time] + 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)
lines(mts[1,, 1:time] - 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)


#Estimando Y e theta

#Amotrando thetha

amostra_theta = array(NA, dim = c(1000, n, time))

#no tempo t = 1
amostra_theta[,,1] = rmvt(1000, delta = mts[,,1], sigma = St[1]/S0*Cts[,,1], df = nt[1])

#no tempo t = 2 em diante
for(i in 2:time){
  if(isSymmetric(Cts[,,i]) == F){
    upper_tri <- upper.tri(Cts[,,i])
    # Tornando a matriz simétrica
    Cts[,,i] <- Cts[,,i] * upper_tri + t(Cts[,,i]) * (1 - upper_tri)
  }
  
  amostra_theta[,,i] = rmvt(1000, delta = mts[,,i], sigma = St[i]/St[i-1]*Cts[,,i], df = nt[i])
}

# retomando as escalas

for(i in 1:time){
  amostra_theta[,,i] = t(t(amostra_theta[,,i])/Ft_padronizado)
}

#amostra do vetor eta

amostra_eta = amostra_theta[,c(2:n),time]

plot(apply(amostra_eta, FUN = quantile, MARGIN = 2, probs = 0.5))
points(eta_ver, col = "red")

apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.025)
apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.5)
apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.975)

#amostra do vetor beta (nas análises seguintes eu pego sempre o último vetor estimado, pois como ele é constante, eu devo pegar o mais atualizado)

#matriz auxiliar
m = matrix(1, nrow = q,
           ncol = n-1)

m[1,] = c(1,0,0)

for (i in 2:(n-1)){
  for (j in 1:(q-1)){
    
    m[j+1,i] = j^(i-1)
    
  }
}

amostra_beta = amostra_eta%*%t(m)


#y estimado
mu = phi = c()

for (i in 1:time){
  mu[i] = t(Ft[,i])%*%mts[,,i]
  phi[i] = qgamma(0.5, nt[i]/2, nt[i]*St[i]/2)
}

y_estimado = qnorm(0.5,mu,sqrt(1/phi))
ic_inf_y = qnorm(0.025,mu,sqrt(1/phi))
ic_sup_y = qnorm(0.975,mu,sqrt(1/phi))

############### Intervalo de Credibilidade de 95% ##########

#Pego o beta mais atualizado, pois ele é constante ao longo do tempo

ic_inf_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.025)
beta_estimado = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.5)
ic_sup_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.975)

#### Y estimado

ggplot(NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y, color = "Dados"), size = 2) +
  geom_line(aes(y = y_estimado, color = "Estimativas"), size = 2) +
  geom_ribbon(aes(ymin = ic_inf_y, ymax = ic_sup_y,
                  fill = "IC 95%"),
              alpha = 0.3) +
  geom_hline(yintercept = 0) +
  labs(x = "Tempo",
       y = "y",
       colour = "") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 45),
        axis.text = element_text(size = 45,
                                 angle = 0),
        legend.text=element_text(size= 30),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_text(size = 45)) +
  scale_color_manual(values = c("Dados" = "black", "Estimativas" = "red"),
                     breaks = c("Estimativas", "Dados")) +
  scale_fill_manual(values = c("IC 95%" = "blue")) +
  scale_y_continuous(limits = c(0.9*min(Y), 1.1*max(Y)))

############ Gráfico dos betas estimados #########
ggplot(NULL,
       aes(x = 0:(q-1))) +
  geom_point(aes(y = beta_estimado), size = 8) +
  geom_errorbar(aes(ymin = ic_inf_beta, ymax = ic_sup_beta),
                width = 0.5, size = 2) +
  geom_point(aes(y = beta_ver, colour = "red"), size = 8, show.legend = F) +
  labs(x = "Lags",
       y = bquote(beta),
       colour = "") +
  theme_hc(base_size = 1) +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        legend.position = NULL,
        panel.grid = element_line(linewidth = 1)) +
  scale_x_continuous(breaks = seq(0, (q-1), by = 1)) +
  scale_y_continuous(breaks = pretty(range(ic_inf_beta, ic_sup_beta)),
                     labels = label_number(decimal.mark = ","))

##### Métricas

mape = sum(abs((Y-y_estimado)/Y))/time * 100
mape

#### Interval Score
interval_score = (ic_sup_y - ic_inf_y) + 2/0.05*(ic_inf_y - Y)*(Y < ic_inf_y) + 
  + 2/0.05*(Y - ic_sup_y)*(Y > ic_sup_y)
sum(interval_score)

######## log_vero

SOMA <- function(MAT,Y,mu){
  soma <- sum(dnorm(Y,mu,sqrt(MAT),log=TRUE))
  return(soma)
}

SOMA(Qt,Y,ft)

##### Curva 2, d = 3 ##########

#contruindo o St,i

#grau do polinômio d = 3

v = c(1:(q-1)) #vetor que multiplica cada matriz (vide as contas na folha)

St0 = apply(X,1,sum)
St1 = X[,-1]%*%v^1
St2 = X[,-1]%*%v^2
St3 = X[,-1]%*%v^3

############## Modelo Linear Dinâmico com Defasagem Polinomial #########

time = nrow(X) #tamanho da série temporal
n = 5        # Dimensão do vetor dos parâmetros a cada tempo
d = 3
 
#vetor F

Ft = matrix(rep(1),n,time)

Ft[2,] = St0
Ft[3,] = St1
Ft[4,] = St2
Ft[5,] = St3

#Bloco da G para o nível

G1 = diag(rep(1,1))

#Bloco da G para as regressoras

G2 = diag(rep(1,n-1)) #n-1 pois tiramos o nível

#matriz G

Gt = bdiag(G1, G2)
Gt = array(as.matrix(Gt), dim = c(n,n,time))
Gt

##### Estimando a série

# Padronizando as variáveis

Ft_padronizado = apply(Ft,1,sd)
Ft_padronizado[1] = 1

Ft = Ft/Ft_padronizado

#Criando os componentes

at = array(NA, dim=c(n,1,time))
Rt = array(NA, dim=c(n,n,time))

mt = array(NA, dim=c(n,1,time))
Ct = array(NA, dim=c(n,n,time))

ft = rep(NA, time)
Qt = rep(NA, time)

et = rep(NA, time)
At = array(NA, dim=c(n,1,time))

St = rep(NA,time)
nt = rep(NA,time)

##### Fator de Desconto

#Bloco do nível

d1 = 0.99
fd1 = 1/(d1)

#Bloco da regressora

d2 = 1
fd2 = 1/(d2)

#Juntando tudo

D_inv = diag(c(fd1, rep(fd2, n-1)))
D_inv

####### Filtro de Kalman

#Passo t  = 0

m0 = c(0, rep(0, n-1))
C0 = diag(c(100, rep(100,n-1)), n, n)

n0 = 3
S0 = 0.5

curve(dgamma(x,n0/2,n0*S0/2), from =0 , to = 20)
abline(v = 1/V_raiz, col = "red")

#Passo t = 1

at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%t(Gt[,,1])) + (D_inv-diag(1, n))*(Gt[,,1]%*%C0%*%t(Gt[,,1]))

ft[1] = t(Ft[,1])%*%at[,,1]
Qt[1] = t(Ft[,1])%*%Rt[,,1]%*%Ft[,1] + S0

et[1] = Y[1] - ft[1]
At[,,1] = Rt[,,1]%*%Ft[,1]/Qt[1]

nt[1] = n0 + 1
St[1] = S0 + S0/nt[1]*(et[1]^2/Qt[1] - 1)

mt[,,1] = at[,,1] + At[,,1]*et[1]
Ct[,,1] = St[1]/S0*(Rt[,,1] - At[,,1]%*%t(At[,,1])*Qt[1])

#Passo t de 2 a time

for (t in 2:time){
  
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] = Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]) + (D_inv-diag(1, n))*(Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]))
  
  ft[t] = t(Ft[,t])%*%at[,,t]
  Qt[t] = t(Ft[,t])%*%Rt[,,t]%*%Ft[,t] + St[t-1]
  
  et[t] = Y[t] - ft[t]
  At[,,t] = Rt[,,t]%*%Ft[,t]/Qt[t]
  
  nt[t] = nt[t-1] + 1
  St[t] = St[t-1] + St[t-1]/nt[t]*(et[t]^2/Qt[t] - 1)
  
  mt[,,t] = at[,,t] + At[,,t]*et[t]
  Ct[,,t] = St[t]/St[t-1]*(Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[t])
  
}

#Analisando

plot.ts(Y)
lines(ft, col = "red")

plot.ts(alpha)
lines(mt[1,,], col = "red")

plot(mt[c(2:n),,time]/Ft_padronizado[-1], col = "red")
points(eta_ver, col = "black")


##### Suavização

#definindo Bt:
Bt = array(NA, dim=c(n,n,time))

mts = mt 
Cts = Ct

for (k in 1:(time-1)){
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time-k+1])
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mts[,,time-k+1] - at[,,time-k+1])
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Cts[,,time-k+1] - Rt[,,time-k+1])%*%t(Bt[,,time-k])
  
  if(is.positive.definite(Cts[,,time-k]) == F){Cts[,,time-k] = make.positive.definite(Cts[,,time-k])}
}

plot.ts(alpha)
lines(mt[1,,1:time], col = "blue")
lines(mt[1,, 1:time] + 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mt[1,, 1:time] - 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mts[1,, 1:time], col = "red")
lines(mts[1,, 1:time] + 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)
lines(mts[1,, 1:time] - 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)


#Amotrando thetha

amostra_theta = array(NA, dim = c(1000, n, time))

#no tempo t = 1
amostra_theta[,,1] = rmvt(1000, delta = mts[,,1], sigma = St[1]/S0*Cts[,,1], df = nt[1])

#no tempo t = 2 em diante
for(i in 2:time){
  if(isSymmetric(Cts[,,i]) == F){
    upper_tri <- upper.tri(Cts[,,i])
    # Tornando a matriz simétrica
    Cts[,,i] <- Cts[,,i] * upper_tri + t(Cts[,,i]) * (1 - upper_tri)
  }
  
  amostra_theta[,,i] = rmvt(1000, delta = mts[,,i], sigma = St[i]/St[i-1]*Cts[,,i], df = nt[i])
}

# retomando as escalas

for(i in 1:time){
  amostra_theta[,,i] = t(t(amostra_theta[,,i])/Ft_padronizado)
}

#amostra do vetor eta

amostra_eta = amostra_theta[,c(2:n),time]

plot(apply(amostra_eta, FUN = quantile, MARGIN = 2, probs = 0.5))
points(eta_ver, col = "red")

apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.025)
apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.5)
apply(amostra_eta, MARGIN = 2, FUN = quantile, probs = 0.975)

#amostra do vetor beta (nas análises seguintes eu pego sempre o último vetor estimado, pois como ele é constante, eu devo pegar o mais atualizado)

#matriz auxiliar
m = matrix(1, nrow = q,
           ncol = n-1)

m[1,] = c(1,0,0,0)

for (i in 2:(n-1)){
  for (j in 1:(q-1)){
    
    m[j+1,i] = j^(i-1)
    
  }
}

amostra_beta = amostra_eta%*%t(m)

apply(amostra_eta, 2, median)

#y estimado
mu = phi = c()

for (i in 1:time){
  mu[i] = t(Ft[,i])%*%mts[,,i]
  phi[i] = qgamma(0.5, nt[i]/2, nt[i]*St[i]/2)
}

y_estimado = qnorm(0.5,mu,sqrt(1/phi))
ic_inf_y = qnorm(0.025,mu,sqrt(1/phi))
ic_sup_y = qnorm(0.975,mu,sqrt(1/phi))

############### Intervalo de Credibilidade de 95% ##########

#Pego o beta mais atualizado, pois ele é constante ao longo do tempo

ic_inf_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.025)
beta_estimado = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.5)
ic_sup_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.975)

#### Y estimado

ggplot(NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y, color = "Dados"), size = 2) +
  geom_line(aes(y = y_estimado, color = "Estimativas"), size = 2) +
  geom_ribbon(aes(ymin = ic_inf_y, ymax = ic_sup_y,
                  fill = "IC 95%"),
              alpha = 0.3) +
  geom_hline(yintercept = 0) +
  labs(x = "Tempo",
       y = "y",
       colour = "") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 45),
        axis.text = element_text(size = 45,
                                 angle = 0),
        legend.text=element_text(size= 30),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_text(size = 45)) +
  scale_color_manual(values = c("Dados" = "black", "Estimativas" = "red"),
                     breaks = c("Estimativas", "Dados")) +
  scale_fill_manual(values = c("IC 95%" = "blue")) +
  scale_y_continuous(limits = c(0.9*min(Y), 1.1*max(Y)))


############ Gráfico dos betas estimados #########
ggplot(NULL,
       aes(x = 0:(q-1))) +
  geom_point(aes(y = beta_estimado), size = 8) +
  geom_errorbar(aes(ymin = ic_inf_beta, ymax = ic_sup_beta),
                width = 0.5, size = 2) +
  geom_point(aes(y = beta_ver, colour = "red"), size = 8, show.legend = F) +
  labs(x = "Lags",
       y = bquote(beta),
       colour = "") +
  theme_hc(base_size = 1) +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        legend.position = NULL,
        panel.grid = element_line(linewidth = 1)) +
  scale_x_continuous(breaks = seq(0, (q-1), by = 1)) +
  scale_y_continuous(breaks = pretty(range(ic_inf_beta, ic_sup_beta)),
                     labels = label_number(decimal.mark = ","))

# log-verossimilhança

vero = SOMA(Qt, Y, y_estimado)
vero

mean(abs((y_estimado-Y)/Y)) * 100 #mape

#### Interval Score
interval_score = (ic_sup_y - ic_inf_y) + 2/0.05*(ic_inf_y - Y)*(Y < ic_inf_y) + 
  + 2/0.05*(Y - ic_sup_y)*(Y > ic_sup_y)
sum(interval_score)


######################## Curva 3 ####################

#vetor de etas verdadeiros:

eta_ver = c(0.1, 2.01, -0.32, 0.012)/50

M_aux = matrix(nrow = q-2,
               ncol = 4)

for (i in 2:(q-1)){
  
  M_aux[i-1,] =  c(1, i, i^2, i^3)
  
}

beta_ver = c(eta_ver[1],
             sum(eta_ver))

for (i in 1:(q-2)){
  beta_ver = c(beta_ver,eta_ver%*%M_aux[i,])
}

plot(x = 1:q, y = beta_ver)
abline(h= 0)

#contruindo o St,i

#grau do polinômio d = 3

v = c(1:(q-1)) #vetor que multiplica cada matriz (vide as contas na folha)

St0 = apply(X,1,sum)
St1 = X[,-1]%*%v^1
St2 = X[,-1]%*%v^2
St3 = X[,-1]%*%v^3

############## Modelo Linear Dinâmico com Defasagem Polinomial #########

time = nrow(X) #tamanho da série temporal
n = 5        # Dimensão do vetor dos parâmetros a cada tempo

#vetor F

Ft = matrix(rep(1),n,time)

Ft[2,] = St0
Ft[3,] = St1
Ft[4,] = St2
Ft[5,] = St3

#Bloco da G para o nível

G1 = diag(rep(1,1))

#Bloco da G para as regressoras

G2 = diag(rep(1,n-1)) #n-1 pois tiramos o nível

#matriz G

Gt = bdiag(G1, G2)
Gt = array(as.matrix(Gt), dim = c(n,n,time))
Gt

#### gerando o Y

V_raiz = sqrt(5) #desvio padrão observacional
W = diag(c(1, rep(0, n-1))) #matrix de variâncial de evolução no tempo 1

alpha_1 = 20

theta_t = matrix(c(alpha_1,eta_ver), nrow = n, ncol = time)

Y = t(Ft[,1])%*%theta_t[,1] + rnorm(1, 0, V_raiz)

for (i in 2:time){
  
  #equação de estados
  theta_t[,i] = Gt[,,i]%*%theta_t[,i-1] + mvrnorm(1, rep(0, n), W)
  
  #equação de evolução
  
  Y[i] = t(Ft[,i])%*%theta_t[,i] + rnorm(1, 0, V_raiz)
  
}

alpha = theta_t[1,]
plot.ts(theta_t[1,])
abline(h=0)

plot.ts(Y)


##### Estimando a série

# Padronizando as variáveis

Ft_padronizado = apply(Ft,1,sd)
Ft_padronizado[1] = 1

Ft = Ft/Ft_padronizado

#Criando os componentes

at = array(NA, dim=c(n,1,time))
Rt = array(NA, dim=c(n,n,time))

mt = array(NA, dim=c(n,1,time))
Ct = array(NA, dim=c(n,n,time))

ft = rep(NA, time)
Qt = rep(NA, time)

et = rep(NA, time)
At = array(NA, dim=c(n,1,time))

St = rep(NA,time)
nt = rep(NA,time)

##### Fator de Desconto

#Bloco do nível

d1 = 0.95
fd1 = 1/(d1)

#Bloco da regressora

d2 = 1
fd2 = 1/(d2)

#Juntando tudo

D_inv = diag(c(fd1, rep(fd2, n-1)))
D_inv

####### Filtro de Kalman

#Passo t  = 0

m0 = c(0, rep(0, n-1))
C0 = diag(c(100, rep(100,n-1)), n, n)

n0 = 3
S0 = 0.5

curve(dgamma(x,n0/2,n0*S0/2), from =0 , to = 20)
abline(v = 1/V_raiz, col = "red")

#Passo t = 1

at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%t(Gt[,,1])) + (D_inv-diag(1, n))*(Gt[,,1]%*%C0%*%t(Gt[,,1]))

ft[1] = t(Ft[,1])%*%at[,,1]
Qt[1] = t(Ft[,1])%*%Rt[,,1]%*%Ft[,1] + S0

et[1] = Y[1] - ft[1]
At[,,1] = Rt[,,1]%*%Ft[,1]/Qt[1]

nt[1] = n0 + 1
St[1] = S0 + S0/nt[1]*(et[1]^2/Qt[1] - 1)

mt[,,1] = at[,,1] + At[,,1]*et[1]
Ct[,,1] = St[1]/S0*(Rt[,,1] - At[,,1]%*%t(At[,,1])*Qt[1])

#Passo t de 2 a time

for (t in 2:time){
  
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] =Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]) + (D_inv-diag(1, n))*(Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t]))
  
  ft[t] = t(Ft[,t])%*%at[,,t]
  Qt[t] = t(Ft[,t])%*%Rt[,,t]%*%Ft[,t] + St[t-1]
  
  et[t] = Y[t] - ft[t]
  At[,,t] = Rt[,,t]%*%Ft[,t]/Qt[t]
  
  nt[t] = nt[t-1] + 1
  St[t] = St[t-1] + St[t-1]/nt[t]*(et[t]^2/Qt[t] - 1)
  
  mt[,,t] = at[,,t] + At[,,t]*et[t]
  Ct[,,t] = St[t]/St[t-1]*(Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[t])
}

#Analisando

plot.ts(Y)
lines(ft, col = "red")

plot.ts(alpha)
lines(mt[1,,], col = "red")

plot(mt[c(2:n),,time]/Ft_padronizado[-1], col = "red")
points(eta_ver, col = "black")


##### Suavização
#definindo Bt:
Bt = array(NA, dim=c(n,n,time))

mts = mt 
Cts = Ct

for (k in 1:(time-1)){
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time-k+1])
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mts[,,time-k+1] - at[,,time-k+1])
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Cts[,,time-k+1] - Rt[,,time-k+1])%*%t(Bt[,,time-k])
  
  if(is.positive.definite(Cts[,,time-k]) == F){Cts[,,time-k] = make.positive.definite(Cts[,,time-k])}
}

plot.ts(alpha)
lines(mt[1,,1:time], col = "blue")
lines(mt[1,, 1:time] + 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mt[1,, 1:time] - 2*sqrt(Ct[1,1,1:time]), col = "blue", lty = 2)
lines(mts[1,, 1:time], col = "red")
lines(mts[1,, 1:time] + 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)
lines(mts[1,, 1:time] - 2*sqrt(Cts[1,1,1:time]), col = "red", lty = 2)


#Amotrando thetha

amostra_theta = array(NA, dim = c(1000, n, time))

#no tempo t = 1
amostra_theta[,,1] = rmvt(1000, delta = mts[,,1], sigma = St[1]/S0*Cts[,,1], df = nt[1])

#no tempo t = 2 em diante
for(i in 2:time){
  if(isSymmetric(Cts[,,i]) == F){
    upper_tri <- upper.tri(Cts[,,i])
    # Tornando a matriz simétrica
    Cts[,,i] <- Cts[,,i] * upper_tri + t(Cts[,,i]) * (1 - upper_tri)
  }
  
  amostra_theta[,,i] = rmvt(1000, delta = mts[,,i], sigma = St[i]/St[i-1]*Cts[,,i], df = nt[i])
}

# retomando as escalas

for(i in 1:time){
  amostra_theta[,,i] = t(t(amostra_theta[,,i])/Ft_padronizado)
}

#amostra do intercepto (alpha) estimado

amostra_alpha = matrix(NA, nrow = 1000, ncol = time)
for (i in 1:time){
  amostra_alpha[,i] = amostra_theta[,1,i]
}

#amostra do vetor eta

amostra_eta = amostra_theta[,c(2:n),time]

plot(apply(amostra_eta, FUN = quantile, MARGIN = 2, probs = 0.5))
points(eta_ver, col = "red")

#amostra do vetor beta (nas análises seguintes eu pego sempre o último vetor estimado, pois como ele é constante, eu devo pegar o mais atualizado)

#matriz auxiliar
m = matrix(1, nrow = q,
           ncol = n-1)

m[1,] = c(1,0,0,0)

for (i in 2:(n-1)){
  for (j in 1:(q-1)){
    
    m[j+1,i] = j^(i-1)
    
  }
}

amostra_beta = amostra_eta%*%t(m)

#y estimado
mu = phi = c()

for (i in 1:time){
  mu[i] = t(Ft[,i])%*%mts[,,i]
  phi[i] = qgamma(0.5, nt[i]/2, nt[i]*St[i]/2)
}

y_estimado = qnorm(0.5,mu,sqrt(1/phi))
ic_inf_y = qnorm(0.025,mu,sqrt(1/phi))
ic_sup_y = qnorm(0.975,mu,sqrt(1/phi))

############### Intervalo de Credibilidade de 95% ##########

#Pego o beta mais atualizado, pois ele é constante ao longo do tempo

ic_inf_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.025)
beta_estimado = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.5)
ic_sup_beta = apply(amostra_beta, MARGIN = 2, FUN = quantile, probs = 0.975)

#### Y estimado

ggplot(NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y, color = "Dados"), size = 2) +
  geom_line(aes(y = y_estimado, color = "Estimativas"), size = 2) +
  geom_ribbon(aes(ymin = ic_inf_y, ymax = ic_sup_y,
                  fill = "IC 95%"),
              alpha = 0.3) +
  geom_hline(yintercept = 0) +
  labs(x = "Tempo",
       y = "y",
       colour = "") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 45),
        axis.text = element_text(size = 45,
                                 angle = 0),
        legend.text=element_text(size= 30),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.y = element_text(size = 45)) +
  scale_color_manual(values = c("Dados" = "black", "Estimativas" = "red"),
                     breaks = c("Estimativas", "Dados")) +
  scale_fill_manual(values = c("IC 95%" = "blue")) +
  scale_y_continuous(limits = c(0.9*min(Y), 1.1*max(Y)))


############ Gráfico dos betas estimados #########
ggplot(NULL,
       aes(x = 0:(q-1))) +
  geom_point(aes(y = beta_estimado), size = 8) +
  geom_errorbar(aes(ymin = ic_inf_beta, ymax = ic_sup_beta),
                width = 0.5, size = 2) +
  geom_point(aes(y = beta_ver, colour = "red"), size = 8, show.legend = F) +
  labs(x = "Lags",
       y = bquote(beta),
       colour = "") +
  theme_hc(base_size = 1) +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        legend.position = NULL,
        panel.grid = element_line(linewidth = 1)) +
  scale_x_continuous(breaks = seq(0, (q-1), by = 1)) +
  scale_y_continuous(breaks = pretty(range(ic_inf_beta, ic_sup_beta)),
                     labels = label_number(decimal.mark = ","))

### Efeito da suavização no nível (não colocado diretamente no TCC)

ggplot(NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y), size = 2) +
  geom_line(aes(y = ft), size = 2, col = "red") +
  geom_line(aes(y = y_estimado), size = 2, col = "blue") +
  geom_hline(yintercept = 0) +
  labs(x = "Tempo",
       y = expression(alpha),
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        legend.position = NULL) +
  scale_y_continuous(limits = c(0.5*min(Y), 1.5*max(Y)))

