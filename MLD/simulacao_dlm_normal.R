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

########### Parte não inclusa no TCC ###########
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