#Pacotes utilizados
setwd("C:/Users/leona/Desktop/Mestrado/1� semestre/Series temporais/Trabalho Computacional")
library(tseries)
library(forecast)
library(seastests)
library(stats)
library(lmtest)
library(Metrics)

#Apagando todas as vari�veis existentes
rm(list = ls(all=TRUE))

#Lendo o Banco de Dados
Consumo <- read.delim2("C:/.../data.txt", header=FALSE)

#Construindo as s�ries 
consumo_ts<- ts(Consumo, frequency=12, start=c(2004,1))

#Checagem sazonalidade
isSeasonal(consumo_ts, test = "wo", freq=12) #WO-test
isSeasonal(consumo_ts, test = "qs", freq=12) #QS test
isSeasonal(consumo_ts, test = "fried", freq=12) #Friedman test
isSeasonal(consumo_ts, test = "kw", freq=12) #Kruskall-Wallis test
isSeasonal(consumo_ts, test = "seasdum", freq=12) #F-test
isSeasonal(consumo_ts, test = "welch", freq=12) #Welch test

#Plotando o gr�fico da S�rie
ts.plot(consumo_ts)
getOption("scipen") 
opt <- options("scipen" = 20) 
ts.plot(consumo_ts, main="Consumo de energia el�trica na rede", xlab="Tempo", ylab="Consumo em MWh") 
options(opt)

#Decomposi��o
consumo_ts_c<- decompose(consumo_ts)
plot(consumo_ts_c)
plot(consumo_ts_c$trend)
plot(consumo_ts_c$seasonal)
plot(consumo_ts_c$random)

#Criando Vari�veis
consumo_in_sample<-ts(Consumo, frequency=12, start=c(2004,1),end=c(2019,5))
consumo_out_of_sample<-consumo_ts[186:197]

#M�todos de Previs�o
consumo1_n <-HoltWinters(consumo_in_sample, gamma = FALSE, beta = FALSE) #naive 
consumo1_h <-HoltWinters(consumo_in_sample, gamma = FALSE) #holt 
consumo1_hwa <-HoltWinters(consumo_in_sample, seasonal=c("additive")) #holt-winters aditivo 
consumo1_hwm <-HoltWinters(consumo_in_sample, seasonal=c("multiplicative")) #holt-winters mutiplicativo

#Checagem de ru�do branco
fit <- consumo_ts_c$random
checkresiduals(fit)
Box.test(consumo_ts_c$random, lag = 1, type = "Ljung-Box", fitdf = 0)
shapiro.test(consumo_ts_c$random)

#Previs�o 12 Passos a Frente
consumo1_n_f <- predict(consumo1_n, n.ahead =12) 
consumo1_h_f <- predict(consumo1_h, n.ahead =12) 
consumo1_hwa_f <- predict(consumo1_hwa, n.ahead =12) 
consumo1_hwm_f <- predict(consumo1_hwm, n.ahead =12) 

#Gr�fico previs�o da s�rie
ts.plot(consumo_out_of_sample , consumo1_n_f, consumo1_h_f , consumo1_hwa_f , consumo1_hwm_f, main="Previs�o da s�rie", xlab="Tempo", ylab="Consumo em Mwh",col=rainbow(6), type="o") 
legend("bottom", inset=-.00,c("Observado","Naive","Holt","Holt-Winters Aditivo","Holt-Winters Multiplicativo"), fill=rainbow(6),cex=0.55)

#Criando fun��o MAPE e c�lculo do MAPE p/ todos os modelos - OUT OF SAMPLE
mape <- function(actual,pred){
  mape <- mean(abs((actual - pred)/actual))*100
  return (mape)
}
#C�lculo do MAPE OUT OF SAMPLE
naive_mape=mape(dados[,3],dados[,4])
holt_mape=mape(dados[,3],dados[,5])
holt_winters_aditivo_mape=mape(dados[,3],dados[,6])
holt_winters_multiplicativo_mape=mape(dados[,3],dados[,7])
#Criando fun��o MAPE e c�lculo do MAPE p/ todos os modelos - IN SAMPLE
MAPE2 <- function(observado,modelo){
  n <- length(observado)
  mape0 <- (sum(abs((residuals(modelo))/observado)))*100
  mape1 <- (1/n)*mape0
  return(mape1)
}
#C�lculo do MAPE In-sample
naive_mape2=MAPE2(consumo_in_sample[1:(length(consumo_in_sample)-1)],consumo1_n)
holt_mape2=MAPE2(consumo_in_sample[1:(length(consumo_in_sample)-2)],consumo1_h)
holt_winters_aditivo_mape2=MAPE2(consumo_in_sample[1:(length(consumo_in_sample)-12)],consumo1_hwa)
holt_winters_multiplicativo_mape2=MAPE2(consumo_in_sample[1:(length(consumo_in_sample)-12)],consumo1_hwm)

#------------------
#BOX AND JENKINS   
#------------------

#Grafico com ACF e PACF
ggtsdisplay(consumo_ts)

lambda <-BoxCox.lambda(consumo_ts)
consumo_ts_bc <- BoxCox(consumo_ts,lambda=lambda)
ggtsdisplay(consumo_ts_bc)

#Diferencia��o n�o sazonal
ndiffs(consumo_ts_bc)
consumo_ts_bc_diff <- diff(consumo_ts_bc,1)
ggtsdisplay(consumo_ts_bc_diff)

#Diferencia��o sazonal
nsdiffs(consumo_ts_bc)
consumo_ts_bc_diff2 <- diff(consumo_ts_bc,12)
ggtsdisplay(consumo_ts_bc_diff2)

#usando o teste de Dickey e Fuller
library(urca)
summary(ur.df(consumo_ts_bc_diff, type='none', lags=0))
summary(ur.df(consumo_ts_bc_diff2, type='none', lags=0))

#Ajustando Arima
dim(consumo_ts)
consumo_ts_treino <- ts(consumo_ts[1:185])#2004~fev.2019
consumo_ts_validacao <- ts(consumo_ts[186:197])#fev.2019~fev.2020

fit <- Arima(y=consumo_ts_treino, order=c(2,1,2), seasonal=list(order=c(1,1,0), 
      period = 12), lambda=TRUE)
#fit <- auto.arima(y=consumo_ts_treino, seasonal=TRUE, lambda=TRUE,stationary = FALSE)
summary(fit)
coeftest(fit)
ts.plot(consumo_ts_treino, xlab="Tempo", ylab="Consumo em Mwh")
lines(fit$fitted, col='blue')

BJ_in_sample=mape(consumo_ts_treino,fit$fitted)

#Realizando forecast para os 12 meses separados para validacao
predi <- forecast(fit,h=12)
ts.plot(predi$mean, main="Previs�o da s�rie", xlab="Tempo", ylab="Consumo em Mwh")

ts.plot(as.numeric(consumo_ts_validacao), type='o',ylim=c(5*10^6, 9*10^6))
lines(as.numeric(predi$mean), col='blue', type="o")

BJ_out_of_sample=mape(as.numeric(consumo_ts_validacao),as.numeric(predi$mean))

#Diagnostico dos res�duos
#tsdiag(fit)

#qqnorm(fit$residuals)
#qqline(fit$residuals)

#Combina��o das previs�es
var_holt_mult=var(residuals(consumo1_hwm))
var_sarima=var(fit$residuals)
k=var_sarima/(var_holt_mult+var_sarima)
previsao<-ts(predi$mean, frequency=12, start=c(2019,6),end=c(2020,5))
previsao_comb=consumo1_hwm_f*(1-k)+previsao*k

#Plotando os modelos de previs�o
observado<-ts(consumo_ts[186:197], frequency=12, start=c(2019,6),end=c(2020,5))
plot(previsao_comb, col="blue", type='o',ylim=c(5*10^6, 9*10^6), ylab= "Consumo de energia em MWh")
lines(observado, col='red', type="o")
legend("bottom", inset=-.00,c("Combina��o","Observado"), fill=c("blue", "red"),cex=0.9)
combinado_mape=mape(as.numeric(previsao_comb),as.numeric(observado))

#Organiza��o dos dados e exporta��o
meses = c("junho", "julho", "agosto", "setembro", "outubro", "novembro" , "dezembro", "janeiro", "fevereiro","mar�o", "abril", "maio")
conversao=data.frame(Dados = consumo_out_of_sample, Naive = consumo1_n_f, 
                     Holt = consumo1_h_f, Holt_Winters_Aditivo = consumo1_hwa_f, 
                     Holt_Winters_Multiplicativo = consumo1_hwm_f, SARIMA = previsao, 
                     Combinado = previsao_comb)

dados = data.frame(Ano = c(rep(2019,7),rep(2020,5)), Meses=meses, 
                   Observado = conversao[,1], Naive = conversao[,2], Holt = conversao[,3],
                   Holt_Winters_Aditivo = conversao[,4], 
                   Holt_Winters_Multiplicativo = conversao[,5], SARIMA = conversao[,6], 
                   Combinado = conversao[,7])

write.csv2(dados, "resultados.csv")

#Organiza��o dos resultados das m�tricas e exporta��o
res_metricas=data.frame(metricas=c("MAPE IN SAMPLE","MAPE OUT OF SAMPLE"), 
                        Naive=c(naive_mape2, naive_mape), Holt=c(holt_mape2, holt_mape), 
                        Holt_Winters_Aditivo=c(holt_winters_aditivo_mape2,
                                               holt_winters_aditivo_mape),
                        Holt_Winters_Multiplicativo=c(holt_winters_multiplicativo_mape2,
                                                      holt_winters_multiplicativo_mape), 
                        Box_and_Jenkins=c(BJ_in_sample,BJ_out_of_sample),
                        Combinado=c("-",combinado_mape))
res_metricas
write.csv2(res_metricas, "metricas.csv")
