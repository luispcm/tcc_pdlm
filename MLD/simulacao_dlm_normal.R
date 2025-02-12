#pacotes
library(tidyverse)
library(epiDisplay)
library(mvtnorm)
library(Matrix)

########## gerando os dados: ###########
set.seed(5)

n = 1 #dimensão de theta
time = 200

Ft = matrix(rep(1),n,time)

#Bloco da G para o nível

G1 = diag(rep(1,1))

#matriz G

Gt = array(1, dim = c(1,1,time))
Gt

#### gerando o Y

V_raiz = sqrt(2) #desvio padrão observacional
W = 0.5 #variância de evolução

alpha_1 = rnorm(1, 20, sqrt(4))

# arbr = c(rgamma(6,1, 1))

theta_t = matrix(c(alpha_1), nrow = n, ncol = time)

Y = t(Ft[,1])%*%theta_t[,1] + rnorm(1, 0, V_raiz)

for (i in 2:time){
  
  #equação de estados
  theta_t[,i] = Gt[,,i]%*%theta_t[,i-1] + rnorm(1, 0, sqrt(W))
  
  #equação de evolução
  
  Y[i] = t(Ft[,i])%*%theta_t[,i] + rnorm(1, 0, V_raiz)
  
}


plot.ts(Y)

###############################################################
###############################################################
## ESTIMANDO A SERIE

n = 1          # Dimensão do vetor dos parâmetros a cada tempo (somente nivel)

# DEFININDO O VETOR F
Ft = matrix(1,nrow=n,ncol=time)

# DEFININDO A MATRIZ G

Gt = array(1, dim=c(n,n,time))

######################### USANDO FATOR DE DESCONTO COM Vt CONHECIDO #####################

##### Filtro de Kalman

### Passo t=0 ###
m0 = rep(0, n)
m0
C0 = diag(1000,n,n)
C0

# Variância observacional
Vt = 2    # Vt conhecido e constante

######### Fator de Desconto

# desconto para o bloco da regressora
delta = 0.99
D = 1/delta

# Definindo as dimensões das matrizes e vetores

at = array(0, dim=c(n,1,time))
Rt = array(0, dim=c(n,n,time))

mt = array(0, dim=c(n,1,time))
Ct = array(0, dim=c(n,n,time))

ft = array(0, dim=c(1,1,time))
Qt = array(0, dim=c(1,1,time))

et = array(0, dim=c(1,1,time))
At = array(0, dim=c(n,1,time))

### Passo t=1 ###

# Priori em t=1
at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%(t(Gt[,,1])))*D

# Previsão 1 passo-a-frente
ft[,,1] = t(Ft[,1])%*%at[,,1]
Qt[,,1] = (t(Ft[,1])%*%Rt[,,1]%*%Ft[,1]) + Vt

# Posteriori em t=1
At[,,1] = (Rt[,,1]%*%Ft[,1])*(1/Qt[,,1])
et[,,1] = Y[1] - ft[,,1]
mt[,,1] = at[,,1] + (At[,,1]*et[,,1])
Ct[,,1] = Rt[,,1] - (At[,,1]%*%(t(At[,,1])))*Qt[,,1]


for(t in 2:time){ # Passo 2 até o passo TIME

  # Priori em t
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] = (Gt[,,t]%*%Ct[,,t-1]%*%(t(Gt[,,t])))*D
  
  # Previsão 1 passo-a-frente
  ft[,,t] = t(Ft[,t])%*%at[,,t]
  Qt[,,t] = (t(Ft[,t])%*%Rt[,,t]%*%Ft[,t]) + Vt
  
  # Posteriori em t
  At[,,t] = (Rt[,,t]%*%Ft[,t])*(1/Qt[,,t])
  et[,,t] = Y[t] - ft[,,t]
  mt[,,t] = at[,,t] + (At[,,t]*et[,,t])
  Ct[,,t] = Rt[,,t] - (At[,,t]%*%(t(At[,,t])))*Qt[,,t]
  
}

ft90 = ft #fator de desconto de 90%
ft95 = ft #fator de desconto de 95%
ft99 = ft #fator de desconto de 99%  

#gráfico
ggplot(data = NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y,
                color = "Y"),
            size = 2) +
  geom_line(aes(y = as.vector(ft90),
                color = "ft90"),
            size = 2) +
  geom_line(aes(y = as.vector(ft95),
                color = "ft95"),
            size = 2) +
  geom_line(aes(y = as.vector(ft99),
                color = "ft99"),
            size = 2) +
  theme_minimal() +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 25),
        legend.position = "bottom") +
  labs(x = "",
       y = "Dados",
       color = "") +
  scale_color_manual(values = c("Y" = "black",
                                "ft90" = "red",
                                "ft95" = "blue",
                                "ft99" = "green"),
                     breaks = c("Y", "ft90", "ft95", "ft99"),
                     labels = c("Dados",
                                expression(delta[theta] ~ "= 0,90"),
                                expression(delta[theta] ~ "= 0,95"),
                                expression(delta[theta] ~ "= 0,99")))



##########################################################################################################################################################################################################################################################################

##########COM A VARIÂNCIA OBSERVACIONAL DESCONHECIDA E FATOR DE DESCONTO#################

############ Filtro de Kalman

### Passo t=0 ###
m0 = rep(0, n)
m0
C0 = diag(1000,n,n)
C0
n0 = 2
S0 = 1

curve(dgamma(x,n0/2,S0/2), from =0 , to = 20)
abline(v = 1/var(Y), col = "red")

############# Utilizando fator de desconto

delta = 0.99
D = 1/delta

########### Definindo as dimensões das matrizes e vetores

at = array(0, dim=c(n,1,time))
Rt = array(0, dim=c(n,n,time))

mt = array(0, dim=c(n,1,time))
Ct = array(0, dim=c(n,n,time))

ft = array(0, dim=c(1,1,time))
Qt = array(0, dim=c(1,1,time))

et = array(0, dim=c(1,1,time))
At = array(0, dim=c(n,1,time))

St = rep(0,time)
nt = rep(0,time)

### Passo t=1

# Priori em t=1
at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = (Gt[,,1]%*%C0%*%(t(Gt[,,1])))*D

# Previsão 1 passo-a-frente
ft[,,1] = t(Ft[,1])%*%at[,,1]
Qt[,,1] = (t(Ft[,1])%*%Rt[,,1]%*%Ft[,1]) + S0

# Posteriori em t=1
At[,,1] = (Rt[,,1]%*%Ft[,1])*(1/Qt[,,1])
et[,,1] = Y[1] - ft[,,1]

#Equações de phi
nt[1] = n0 + 1
St[1] = S0 + (S0/nt[1])*(((et[,,1])^2/Qt[,,1])-1)

mt[,,1] = at[,,1] + (At[,,1]*et[,,1])
Ct[,,1] = (St[1]/S0)*(Rt[,,1] - At[,,1]%*%t(At[,,1])*Qt[,,1])


for(t in 2:time){ # Passo 2 até o passo TIME
  
  # Priori em t
  at[,,t] = Gt[,,t]%*%mt[,,t-1]
  Rt[,,t] = (Gt[,,t]%*%Ct[,,t-1]%*%(t(Gt[,,t])))*D
  
  # Previsão 1 passo-a-frente
  ft[,,t] = t(Ft[,t])%*%at[,,t]
  Qt[,,t] = (t(Ft[,t])%*%Rt[,,t]%*%Ft[,t]) + St[t-1]
  
  # Posteriori em t
  At[,,t] = (Rt[,,t]%*%Ft[,t])%*%(1/Qt[,,t])
  et[,,t] = Y[t] - ft[,,t]
  
  #Equações de phi
  nt[t] = nt[t-1] + 1
  St[t] = St[t-1] + ((St[t-1]/nt[t])*(((et[,,t])^2/Qt[,,t])-1))
  
  mt[,,t] = at[,,t] + (At[,,t]*et[,,t])
  Ct[,,t] = St[t]/St[t-1]*(Rt[,,t] - At[,,t]%*%t(At[,,t])*Qt[,,t])
  
}

ft90 = ft #fator de desconto de 90%
ft95 = ft #fator de desconto de 95%
ft99 = ft #fator de desconto de 99%  

mt90 = mt #fator de desconto de 90%
mt95 = mt #fator de desconto de 95%
mt99 = mt #fator de desconto de 99%  

#gráfico
ggplot(data = NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y,
                color = "Y"),
            size = 2) +
  geom_line(aes(y = as.vector(ft90),
                color = "ft90"),
            size = 2) +
  geom_line(aes(y = as.vector(ft95),
                color = "ft95"),
            size = 2) +
  geom_line(aes(y = as.vector(ft99),
                color = "ft99"),
            size = 2) +
  theme_minimal() +
  theme(axis.title = element_text(size = 45),
        axis.text = element_text(size = 45),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 25),
        legend.position = "bottom") +
  labs(x = "",
       y = "Dados",
       color = "") +
  scale_color_manual(values = c("Y" = "black",
                                "ft90" = "red",
                                "ft95" = "blue",
                                "ft99" = "green"),
                     breaks = c("Y", "ft90", "ft95", "ft99"),
                     labels = c("Dados",
                                expression(delta[theta] ~ "= 0,90"),
                                expression(delta[theta] ~ "= 0,95"),
                                expression(delta[theta] ~ "= 0,99")))

# gráfico do theta
ggplot(data = NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = theta_t[1,],
                color = "Y"),
            size = 2) +
  geom_line(aes(y = as.vector(mt90[1,,]),
                color = "mts90"),
            size = 2) +
  geom_line(aes(y = as.vector(mt95[1,,]),
                color = "mts95"),
            size = 2) +
  geom_line(aes(y = as.vector(mt99[1,,]),
                color = "mts99"),
            size = 2) +
  theme_minimal() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text=element_text(size= 25),
        legend.title = element_text(size = 20),
        legend.position = "bottom") +
  labs(x = "",
       y = expression(theta),
       color = "") +
  scale_color_manual(values = c("Y" = "black",
                                "mts90" = "red",
                                "mts95" = "blue",
                                "mts99" = "green"),
                     labels = c(expression(theta),
                                "Fator de Desconto de 90%",
                                "Fator de Desconto de 95%",
                                "Fator de Desconto de 99%"))

###############################################################################################################################################################################################

########### NÃO VOU COLOCAR ESSA PARTE NO TCC ###########
########### UTILIZANDO A SUAVIZAÇÃO ###############

#definindo Bt:
Bt = array(0, dim=c(n,n,time))

mts = mt #isso porque para k = 0, o mts[time] = mt[time], ou seja, no último tempo, mts e mt sãos iguais
Cts = Ct #o mesmo comentário de cima é válido para esse

for (k in 1:(time-1)){
  
  # k equivalente time-k         #Gt[,,time-k+1]
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time+1-k])
  #no livro, at(-k) = mt no tempo t-k e Rt(-k) é Ct no mesmo tempo citado
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mt[,,time-k+1] - at[,,time-k+1])
  
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Ct[,,time-k+1]-Rt[,,time-k+1])%*%t(Bt[,,time-k])
}

mts90 = mts #fator de desconto de 90%
mts95 = mts #fator de desconto de 95%
mts99 = mts #fator de desconto de 99%  

#gráfico
ggplot(data = NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = theta_t[1,],
                color = "theta"),
            size = 2) +
  geom_line(aes(y = as.vector(mts90[1,,]),
                color = "mts90"),
            size = 2) +
  geom_line(aes(y = as.vector(mts95[1,,]),
                color = "mts95"),
            size = 2) +
  geom_line(aes(y = as.vector(mts99[1,,]),
                color = "mts99"),
            size = 2) +
  theme_minimal() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        legend.text=element_text(size= 20),
        legend.title = element_text(size = 25),
        legend.position = "bottom") +
  labs(x = "",
       y = expression(theta),
       color = "") +
  scale_color_manual(values = c("theta" = "black",
                                "mts90" = "red",
                                "mts95" = "blue",
                                "mts99" = "green"),
                     labels = c(expression(theta),
                                "Fator de Desconto de 90%",
                                "Fator de Desconto de 95%",
                                "Fator de Desconto de 99%"))

#################### Modelo com Suavização e Sazonalidade ###################

n = 7         # Dimensão do vetor dos parâmetros a cada tempo (somente nivel)

Ft = matrix(rep(1),n,time)
Ft[c(2:n),] = c(rep(c(1,0), 3)) #respectivo a sazonalidade

#Bloco da G para o nível

G1 = diag(rep(1,1))

# Bloco da G para a sazonalidade, com h = (p-1)/2 por p ser ímpar

w1=(2*pi)/(7/1)
G3=matrix(c(cos(w1),-sin(w1),sin(w1),cos(w1)),2,2)
G3
w2=(2*pi)/(7/2)
G4=matrix(c(cos(w2),-sin(w2),sin(w2),cos(w2)),2,2)
G4
w3=(2*pi)/(7/3)
G5=matrix(c(cos(w3),-sin(w3),sin(w3),cos(w3)),2,2)
G5

#matriz G

Gt = bdiag(G1, G3, G4, G5)
Gt = array(as.matrix(Gt), dim = c(n,n,time))
Gt

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
fd1 = 1/sqrt(d1)

#Bloco da sazonalidade

d3 = 0.99
fd3 = 1/sqrt(d3)

#Juntando tudo

D = diag(c(fd1, rep(fd3, n-1)))
D

####### Filtro de Kalman

#Passo t  = 0

m0 = rep(0, n)
C0 = diag(1000, n, n)

n0 = 1
S0 = 0.5

curve(dgamma(x,n0/2,n0*S0/2), from =0 , to = 20)
abline(v = 1/V_raiz, col = "red")

#Passo t = 1

at[,,1] = Gt[,,1]%*%m0
Rt[,,1] = D%*%Gt[,,1]%*%C0%*%t(Gt[,,1])%*%t(D)

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
  Rt[,,t] = D%*%Gt[,,t]%*%Ct[,,t-1]%*%t(Gt[,,t])%*%t(D)
  
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

##### Suavização

#definindo Bt:
Bt = array(NA, dim=c(n,n,time))

mts = mt 
Cts = Ct

for (k in 1:(time-1)){
  Bt[,,time-k] = Ct[,,time-k]%*%t(Gt[,,time-k+1])%*%solve(Rt[,,time-k+1], tol = 1e-18)
  mts[,,time-k] = mt[,,time-k] + Bt[,,time-k]%*%(mts[,,time-k+1] - at[,,time-k+1])
  Cts[,,time-k] = Ct[,,time-k] + Bt[,,time-k]%*%(Cts[,,time-k+1] - Rt[,,time-k+1])%*%t(Bt[,,time-k])
}

alpha = theta_t[1,]
plot.ts(alpha)
lines(mts[1,,], col = "red")

y_estimado_90 = c(NA)
y_estimado_95 = c(NA)
y_estimado_99 = c(NA)

for (t in 1:time){
  y_estimado_90[t] = t(Ft[,t])%*%mts[,,t]
}

for (t in 1:time){
  y_estimado_95[t] = t(Ft[,t])%*%mts[,,t]
}

for (t in 1:time){
  y_estimado_99[t] = t(Ft[,t])%*%mts[,,t]
}

#gráfico
ggplot(data = NULL,
       aes(x = 1:time)) +
  geom_line(aes(y = Y,
                color = "theta"),
            size = 2) +
  geom_line(aes(y = y_estimado_90,
                color = "mts90"),
            size = 2) +
  geom_line(aes(y = y_estimado_95,
                color = "mts95"),
            size = 2) +
  geom_line(aes(y = y_estimado_99,
                color = "mts99"),
            size = 2) +
  theme_minimal() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        legend.text=element_text(size= 20),
        legend.title = element_text(size = 25),
        legend.position = "bottom") +
  labs(x = "",
       y = expression(theta),
       color = "") +
  scale_color_manual(values = c("theta" = "black",
                                "mts90" = "red",
                                "mts95" = "blue",
                                "mts99" = "green"),
                     labels = c(expression(theta),
                                "Fator de Desconto de 90%",
                                "Fator de Desconto de 95%",
                                "Fator de Desconto de 99%"))

