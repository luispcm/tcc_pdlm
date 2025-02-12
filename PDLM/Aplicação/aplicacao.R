library(pacman)

p_load("tidyverse",
       "reshape2",
       "microdatasus",
       "mvtnorm")

#########################################################
memory.limit(9999999999)

SP_bruto_2019 = fetch_datasus(year_start = 2018,
                              month_start = 1,   #Tirei o primeiro mês para não pegar ano epidemiológico de 2017 e tirei o mês de outubro a dezembro de 2019 porque há uma queda brusca que pode ser um erro
                              year_end = 2019,
                              month_end = 12,
                              uf = "SP",
                              information_system = "SIH-RD",
                              vars = c("UF_ZI","MUNIC_RES", "DIAG_PRINC", "DIAG_SECUN", "IDADE", "NASC",
                                       "COD_IDADE", "DT_INTER", "DT_SAIDA", "MORTE"))

SP_19_Resp = SP_bruto_2019 |> 
  filter(grepl("^J", DIAG_PRINC), MUNIC_RES == "355030") |> 
  mutate(DT_INTER = as.Date(DT_INTER, format =  "%Y%m%d"),
         ano_epi = epiyear(DT_INTER)) |> 
  filter(ano_epi %in% c(2018,2019)) |> 
  group_by(DT_INTER) |> 
  summarise(casos_resp = n(),
            obitos_resp = sum(MORTE, na.rm = TRUE))

SP_19_Circ = SP_bruto_2019 |> 
  filter(grepl("^I", DIAG_PRINC), MUNIC_RES == "355030") |> 
  mutate(DT_INTER = as.Date(DT_INTER, format =  "%Y%m%d"),
         ano_epi = epiyear(DT_INTER)) |> 
  filter(ano_epi %in% c(2018,2019)) |> 
  group_by(DT_INTER) |> 
  summarise(casos_resp = n(),
            obitos_resp = sum(MORTE, na.rm = TRUE))

par(mfrow=c(1,2))
plot(SP_19_Resp$casos_resp, type="l")
plot(SP_19_Resp$obitos_resp, type="l")


#manipulando a base de dados do SISAM

sisam = rbind(task_9045_dados_sisam_2019, task_9045_dados_sisam_2018)

#filtrando pela cidade de SP e guardando no objeto "base"

base = sisam %>% 
  filter(uf_nome == "SÃO PAULO",
         municipio_nome %in% c("CAJAMAR", "CAIEIRAS", "MAIRIPORÃ", "GUARULHOS", 
                               "ITAQUAQUECETUBA", "POÁ", "FERRAZ DE VASCONCELOS", 
                               "MAUÁ", "SANTO ANDRÉ", "SÃO CAETANO DO SUL", "DIADEMA", 
                               "SÃO BERNARDO DO CAMPO", "SÃO VICENTE", "ITANHAÉM", "COTIA", 
                               "JUQUITIBA", "EMBU-GUAÇU", "ITAPECERICA DA SERRA", "EMBU DAS ARTES", 
                               "TABOÃO DA SERRA", "OSASCO", "BARUERI", "SANTANA DE PARNAÍBA"))

#pegando apenas a data da coleta
base$datahora = substr(base$datahora, 1, 10)

#transformando na classe data
base$datahora = as.Date(base$datahora, format = "%Y-%m-%d")

#filtrando pelo ano epidemiológico de 2019 e resumindo a base por data
base = base %>% 
  mutate(ano_epi = epiyear(datahora)) |> 
  filter(ano_epi %in% c(2018,2019)) |> 
  group_by(datahora) |>
  summarize(co_ppb_medio = mean(co_ppb),
            no2_ppb_medio = mean(no2_ppb),
            o3_ppb_medio = mean(o3_ppb),
            pm25_ugm3_medio = mean(pm25_ugm3, na.rm = T),
            so2_ugm3_medio = mean(so2_ugm3),
            precipitacao_mmdia_media = mean(precipitacao_mmdia),
            temperatura_c_media = mean(temperatura_c),
            umidade_relativa_percentual_media = mean(umidade_relativa_percentual),
            vento_direcao_grau_medio = mean(vento_direcao_grau),
            vento_velocidade_ms_media = mean(vento_velocidade_ms))


# Juntando a base do SUS com a do SISAM

names(base)[1] = "DT_INTER" #renomeando a primeira coluna para fazer o left_join por meio dela

base_circ = left_join(SP_19_Circ, base)
base_resp = left_join(SP_19_Resp, base)

###### Selecionando o período de tempo

base_circ <- base_circ[base_circ$DT_INTER >= as.Date("2019-01-01") & base_circ$DT_INTER <= as.Date("2019-11-30"), ]

base_resp <- base_resp[base_resp$DT_INTER >= as.Date("2019-01-01") & base_resp$DT_INTER <= as.Date("2019-11-30"), ]

# manipulando os dados populacionais da cidade de SP

populacao = apply(populacao_sp[,-1], 2, sum) #população de SP por ano (de 2011 até 2019)

populacao_2018 = populacao[8] #população de SP em 2018
populacao_2019 = populacao[9] #população de SP em 2019

################## calculando a taxa de internações por 100 mil habitantes

# para doenças do aparelho circulatório

base_circ$taxa_inter[substr(base_circ$DT_INTER,1,4) == 2018] = base_circ$casos_resp[substr(base_circ$DT_INTER,1,4) == 2018]/populacao_2018 * 100000

base_circ$taxa_inter[substr(base_circ$DT_INTER,1,4) == 2019] = base_circ$casos_resp[substr(base_circ$DT_INTER,1,4) == 2019]/populacao_2019 * 100000

# para doenças do aparelho respiratório

base_resp$taxa_inter[substr(base_resp$DT_INTER,1,4) == 2018] = base_resp$casos_resp[substr(base_resp$DT_INTER,1,4) == 2018]/populacao_2018 * 100000

base_resp$taxa_inter[substr(base_resp$DT_INTER,1,4) == 2019] = base_resp$casos_resp[substr(base_resp$DT_INTER,1,4) == 2019]/populacao_2019 * 100000

########### calculando a taxa de óbitos por 100 mil habitantes

# para doenças do aparelho circulatório

base_circ$taxa_obitos[substr(base_circ$DT_INTER,1,4) == 2018] = base_circ$obitos_resp[substr(base_circ$DT_INTER,1,4) == 2018]/populacao_2018 * 100000

base_circ$taxa_obitos[substr(base_circ$DT_INTER,1,4) == 2019] = base_circ$obitos_resp[substr(base_circ$DT_INTER,1,4) == 2019]/populacao_2019 * 100000

# para doenças do aparelho respiratório

base_resp$taxa_obitos[substr(base_resp$DT_INTER,1,4) == 2018] = base_resp$obitos_resp[substr(base_resp$DT_INTER,1,4) == 2018]/populacao_2018 * 100000

base_resp$taxa_obitos[substr(base_resp$DT_INTER,1,4) == 2019] = base_resp$obitos_resp[substr(base_resp$DT_INTER,1,4) == 2019]/populacao_2019 * 100000

########### Análise de dados ##########

############### gráfico do pm2.5:

ggplot(base_resp,
       aes(x = DT_INTER)) +
  geom_line(aes(y = pm25_ugm3_medio),
            size = 1) +
  labs(x = "Mês-Ano",
       y = expression(paste("PM2,5 (", mu, "g/m"^3, ")")),
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 38),
        axis.text = element_text(size = 30),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5)) +
  scale_x_date(labels = scales::date_format("%b-%Y"),
               breaks = seq(from = min(base_resp$DT_INTER),
                            to = max(base_resp$DT_INTER),
                            by = "1 month"))


############### gráfico da taxa de internados (casos respiratórios):

ggplot(base_resp,
       aes(x = DT_INTER)) +
  geom_line(aes(y = taxa_inter),
            size = 1) +
  labs(x = "Mês-Ano",
       y = "T. de internações por 100 mil hab.",
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 38),
        axis.text = element_text(size = 30),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.title.y = element_text(hjust = 1,
                                    size = 30)) +
  scale_x_date(labels = scales::date_format("%b-%Y"),
               breaks = seq(from = min(base_resp$DT_INTER),
                            to = max(base_resp$DT_INTER),
                            by = "1 month"))

############### gráfico da umidade:

ggplot(base_resp,
       aes(x = DT_INTER)) +
  geom_line(aes(y = umidade_relativa_percentual_media),
            size = 1) +
  labs(x = "Mês-Ano",
       y = "Umidade Relativa do Ar Média Diária (%)",
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 38),
        axis.text = element_text(size = 30),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.title.y = element_text(hjust = 1,
                                    size = 30)) +
  scale_x_date(labels = scales::date_format("%b-%Y"),
               breaks = seq(from = min(base_resp$DT_INTER),
                            to = max(base_resp$DT_INTER),
                            by = "1 month"))

############### gráfico da temperatura:

ggplot(base_resp,
       aes(x = DT_INTER)) +
  geom_line(aes(y = temperatura_c_media),
            size = 1) +
  labs(x = "Mês-Ano",
       y = "Temperatura Média Diária em °C",
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 38),
        axis.text = element_text(size = 30),
        legend.text=element_text(size= 40),
        legend.title = element_text(size = 45),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.title.y = element_text(hjust = 1,
                                    size = 30)) +
  scale_x_date(labels = scales::date_format("%b-%Y"),
               breaks = seq(from = min(base_resp$DT_INTER),
                            to = max(base_resp$DT_INTER),
                            by = "1 month"))

########## gráfico padronizado

ggplot(base_resp,
       aes(x = DT_INTER)) +
  geom_line(aes(y = pm25_ugm3_medio/max(pm25_ugm3_medio),
                color = "pm25"),
            size = 1) +
  geom_line(aes(y = umidade_relativa_percentual_media/(max(umidade_relativa_percentual_media)),color = "umidade"),
            size = 1) +
  geom_line(aes(y = taxa_inter/(max(taxa_inter)),color = "Y"),
            size = 1) +
  geom_line(aes(y = temperatura_c_media/(max(temperatura_c_media)),color = "temperatura"),
            size = 1) +
  labs(x = "Mês-Ano",
       y = "",
       colour = "") +
  theme_minimal() +
  theme(axis.title = element_text(size = 38),
        axis.text = element_text(size = 30),
        legend.text=element_text(size= 20),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.title.y = element_text(hjust = 1),
        legend.position = "bottom") +
  scale_x_date(labels = scales::date_format("%m-%Y"),
               breaks = seq(from = min(base_resp$DT_INTER),
                            to = max(base_resp$DT_INTER),
                            by = "1 month")) +
  scale_color_manual(values = c("pm25" = "green",
                                "Y" = "black",
                                "umidade" = "darkblue",
                                "temperatura" = "red"),
                     labels = c(expression(paste("PM2,5 (", mu, "g/m"^3,")   ")),
                                "Temperatura (°C)   ",
                                "Umidade Relativa do Ar (%)   ",
                                "Taxa de internação por 100 mil hab."))


####################### DPLDM (Sitema Respiratório - Internações) #####################

memory.limit(9999999999)

########################################## M1
DPLDM(Y = base_resp$taxa_inter, x = base_resp$DT_INTER, regressora = base_resp$pm25_ugm3_medio, lags = 10, d = 2, fd_nivel = 1, fd_phi = 1)

DPLDM(Y = base_resp$taxa_inter, x = base_resp$DT_INTER, regressora = base_resp$pm25_ugm3_medio, lags = 10, d = 3, fd_nivel = 1, fd_phi = 1)

########################################## M2

DPLDM(Y = base_resp$taxa_inter, x = base_resp$DT_INTER, regressora = base_resp$pm25_ugm3_medio, lags = 10, d = 2, fd_nivel = 0.96, fd_phi = 0.99)

DPLDM(Y = base_resp$taxa_inter, x = base_resp$DT_INTER, regressora = base_resp$pm25_ugm3_medio, lags = 10, d = 3, fd_nivel = 0.96, fd_phi = 0.99)

d########################################## M3
### Considerando só umidade
# DPLDM(Y = base_resp$taxa_inter,
#       x = base_resp$DT_INTER,
#       regressora1 = base_resp$pm25_ugm3_medio,
#       regressora2 = base_resp$umidade_relativa_percentual_media,
#       lags = 10,
#       lag_umidade = 7,
#       d = 3,
#       fd_nivel = 0.98,
#       fd_phi = 0.99,
#       fd_umidade = 1,
#       umidade_corte = 90)

### Considerando só temperatura
DPLDM(Y = base_resp$taxa_inter,
      x = base_resp$DT_INTER,
      regressora1 = base_resp$pm25_ugm3_medio,
      regressora2 = base_resp$temperatura_c_media,
      lags = 10,
      lag_temperatura = 0,
      d = 3,
      fd_nivel = 0.96,
      fd_phi = 0.99,
      fd_temperatura = 1,
      quantil_corte = 0.85)

