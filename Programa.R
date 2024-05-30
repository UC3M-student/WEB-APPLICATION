library(fpp2)
library(tseries)
library(quantmod)
library(stats)
library("goftest")
library(plyr)

#### FUNCIONES PROPIEAS ####

diagnosys_phase = function(modelo,time_serie_data) {
  
  problems = 0 #Numero de errores en el diagnostico
  
  # COEFICIENTES NO SIGNIFICATIVOS
  
  coef_no_significativos = 0
  
  tryCatch( {
    
    for (i in 1:length(modelo$coef)) {
      
      if (qnorm(c(0.025,0.975),0,(modelo$var.coef^0.5)[i,i])[1] < modelo$coef[i] &  qnorm(c(0.025,0.975),0,(modelo$var.coef^0.5)[i,i])[2] > modelo$coef[i] ) {
        
        problems = problems + 1 
        
        coef_no_significativos = coef_no_significativos + 1
        
      }
      
    }
    
  }, error = function(e) {
    
    problems = 1000
    return(problems)
    
  }
  
  )
  
  # RESIDUOS INCORRELACIONADOS
  
  residuos_incorrelacionados = 0
  
  res <- residuals(modelo) #Residuals
  ggacf_tsd = ggAcf(time_serie_data)
  
  incorrelacion_pvalue = Box.test(res, lag = max(ggacf_tsd$data$lag), type = c("Ljung-Box"))
  
  residuos_incorrelacionados = incorrelacion_pvalue$p.value
  
  # if (incorrelacion$p.value < 0.05) {
  #   
  #   problems = problems + 1
  #   
  #   residuos_incorrelacionados = residuos_incorrelacionados + 1
  #   
  # }
  
  ### NORMALIDAD EN LOS RESIDUOS
  
  normalidad_residuos = 0
  
  CvM_normal = cvm.test(res,"pnorm",0,sd(res),estimated = TRUE) # H0: Normalidad
  ad_normal = ad.test(res,"pnorm",0,sd(res),estimated = TRUE)
  jb = jarque.bera.test(res)
  
  normalidad_residuos = CvM_normal$p.value
  
  if (CvM_normal$p.value < 0.05) {
    
    problems = problems + 1
    #normalidad_residuos = normalidad_residuos + 1
    
  }
  
  # Media marginal igual a cero - Esperanza nula
  
  media_marginal_cero = 1
  
  if (qnorm(c(0.025,0.975),0,sd(res)/sqrt(length(res)))[1] > mean(res) | qnorm(c(0.025,0.975),0,sd(res)/sqrt(length(res)))[2] < mean(res)) {
    
    problems = problems + 1
    media_marginal_cero = 0
    
  }
  
  
  possible_solutions_table <- matrix(
    c(problems, coef_no_significativos, residuos_incorrelacionados,
      normalidad_residuos, media_marginal_cero),
    ncol = 5
  )
  
  colnames(possible_solutions_table) <- c(
    "Diagnosys_Problems", "Coef_no_sign.",
    "Incorrelacion p.value", "Normalidad p.value", "Esperanza nula"
  )
  
  return(possible_solutions_table)
  
} 

mean_stationarity = function(time_serie_data) {
  
  options(warn=-1) #Suppress warnings
  
  library(tseries)
  
  differences = 0
  
  for (i in 1:30) {
    
    kpss_ = kpss.test(time_serie_data,null = "Level")
    pp_ =  pp.test(time_serie_data)
    adf_ = adf.test(time_serie_data)
    
    if (kpss_$p.value > 0.05 & pp_$p.value < 0.05 | kpss_$p.value < 0.05 & pp_$p.value < 0.05) {
      
      
      options(warn=0) #Undo suppress warnings
      
      return(differences)
      
      #stop()
      
    }
    
    
    time_serie_data = diff(time_serie_data, differences = i)
    
    differences = differences + 1
    
  }
  
} 

ARIMA_Summary_best_fit <- function(time_serie_data) {
  
  count_problems_vector <- c()
  ar_vector <- c()
  difference_vector <- c()
  ma_vector <- c()
  BIC_vector <- c()
  AIC_vector <- c()
  coef_no_significativos <- c()
  residuos_incorrelacionados <- c()
  normalidad_residuos <- c()
  media_marginal_cero <- c()
  var_epsilon = c()
  
  contador <- 1
  
  for (i in 0:4) { # Ar
    for (j in 0:4) { # Ma
      for (k in 0:4) { # difference
        
        if (j == 0 & i == 0) {
          next
        } else {
          tryCatch({
            fit_ <- Arima(time_serie_data, c(i, k, j), include.constant = FALSE)
            count_problems_vector[contador] <- diagnosys_phase(fit_, time_serie_data)[1]
            ar_vector[contador] <- i
            ma_vector[contador] <- j
            difference_vector[contador] <- k
            BIC_vector[contador] <- fit_$bic
            AIC_vector[contador] <- fit_$aic
            coef_no_significativos[contador] <- diagnosys_phase(fit_, time_serie_data)[2]
            residuos_incorrelacionados[contador] <- diagnosys_phase(fit_, time_serie_data)[3]
            normalidad_residuos[contador] <- diagnosys_phase(fit_, time_serie_data)[4]
            media_marginal_cero[contador] <- diagnosys_phase(fit_, time_serie_data)[5]
            var_epsilon[contador] = fit_$sigma2 #error
            
            contador <- contador + 1
          }, error = function(e) {
            # Manejar el error aquí
            #next
            #cat("Error ajustando modelo ARIMA con parámetros:", i, k, j, "\n")
          })
        }
      }
    }
  }
  
  possible_solutions_table <- matrix(
    c(ar_vector, difference_vector, ma_vector, AIC_vector, BIC_vector,
      count_problems_vector, coef_no_significativos, residuos_incorrelacionados,
      normalidad_residuos, media_marginal_cero,var_epsilon),
    ncol = 11
  )
  colnames(possible_solutions_table) <- c(
    "p", "d", "q", "AIC", "BIC", "Diagnosys_Problems", "Coef_no_signi",
    "Incorrelacion", "Normalidad", "Media_marginal_cero", "Var_Epsilon - Sigma^2"
  )

  
  mean_stationarity_output = mean_stationarity(time_serie_data)
  
  df_possible_solutions_table = as.data.frame(possible_solutions_table) 
  
  df_possible_solutions_table = df_possible_solutions_table[with(df_possible_solutions_table, order(Diagnosys_Problems, d,BIC)),]
  
  df_possible_solutions_table = df_possible_solutions_table[grep(mean_stationarity_output, df_possible_solutions_table$d),][1,]
  
  return(df_possible_solutions_table)
  
}

###########################

#### DATOS ####

Mortality_South_Korea = read.csv2("https://raw.githubusercontent.com/UC3M-student/Proyecto_Series_Temporales/main/Mortality_South_Korea.csv", skip = 1)

y_45<-ts(as.numeric(Mortality_South_Korea$X45))
y_46<-ts(as.numeric(Mortality_South_Korea$X46))
y_47<-ts(as.numeric(Mortality_South_Korea$X47))
y_48<-ts(as.numeric(Mortality_South_Korea$X48))
y_49<-ts(as.numeric(Mortality_South_Korea$X49))
y_50<-ts(as.numeric(Mortality_South_Korea$X50))
y_51<-ts(as.numeric(Mortality_South_Korea$X51))
y_52<-ts(as.numeric(Mortality_South_Korea$X52))
y_53<-ts(as.numeric(Mortality_South_Korea$X53))

#######################

#### COMIENZO DEL PROGRAMA ###
my_ARIMA_prediccion = ARIMA_Summary_best_fit(y_45)

modelo_ARIMA = arima(y_45, c(my_ARIMA_prediccion[1,1],my_ARIMA_prediccion[1,2],my_ARIMA_prediccion[1,3]))


for (i in 1:10) {
  
  fcast <- forecast(modelo_ARIMA, h = i)
  one<-rnorm(100000,fcast$mean[i], ((unname(fcast$upper[i,2]) - fcast$mean[i])/ (1.644854) )) #z_score_95 <- 1.644854
  hist(one, freq = F, main= paste("Future moment", i) , col = "lightblue",  xlab = "Value", ylab = "Frequency")
  
}


table_predicter = as.data.frame(fcast)[,c(1,4,5)]

rownames(table_predicter) = c("Future Moment 1","Future Moment 2","Future Moment 3","Future Moment 4","Future Moment 5",
                              "Future Moment 6","Future Moment 7","Future Moment 8","Future Moment 9","Future Moment 10")





