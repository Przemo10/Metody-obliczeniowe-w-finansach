#' Model Hestona 1 iteracja
#' Funkcja \code{Heston.model.one.iteration} wyznacza cenę 
#' @param n liczba trajektorii procesu Wienera
#' @param S0 wartość procesu S0 w chwili 0
#' @param V0 wartość V0
#' @param r stopa procentowa
#' @param sigma zmienność
#' @param a parametr a modelu Hestona
#' @param b parametr b modelu Hestona 
#' @param rho korelacja procesów Wienera
#' @param Time czas
#' @param K cena wykonania
#' 
#' @details Trajektorie procesu V zostały wyznaczone z wykorzystaniem schematu Mildsteina.
#' Trajektorie procesu S z użyciem schematu Eulera z zastąpieniem vt = max(vt,0)(full truncation scheme).
#' 
#' 
#' @export

Heston.model.one.iteration <- function(n, S0 = 100, V0 = 0.09, sigma = 1.0, rho = -0.3, K = 100, b =0.09,  delta.t = 0.005 , alpha1=  0.01, alpha2 = 0.07071068, theta1= 1.00025,df = 0.7788008){
  
  ############################################################################
  ################ GENEROWANIE ZMIENNYCH NORMALNYCH O KORELACJI rho ##########
  ############################################################################
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  ZV <- Z1
  ZS <- rho * Z1 + sqrt(1- rho^2) * Z2
  ##########################################################
  #### GENEROWANIE TRAJEKTORII Vt - SCHEMAT MILDSTEINA #####
  ##########################################################
  
  V <- numeric(n)
  V[1] <- V0
  beta <- sigma^2 * delta.t * (ZV^2 -1 )/4 
  
  ZS2 <- ZS * sqrt(delta.t)
  for ( i in 1:(n-1)){
    
    V[i+1] <- V[i] + alpha1 * (b - max(V[i],0)) + alpha2 * sqrt(max(V[i],0)) * ZV[i] + beta[i]
  }
  
  #####################################################################
  ######## WYZNACZENIE WARTOŚĆI PROCESU S W CHWLI T(schemat Eulera)####
  #####################################################################
  Vprim <- ifelse(V>0,V,0)
  S <- S0 * prod(theta1+ sqrt(Vprim)  * ZS2)
  
  #####################################################################
  ###### WYZNACZENIE WARTOŚCI OPCJI EUROPEJSKIEJ KUPNA W T=0###########
  #####################################################################
  
  df * max(S- K ,0)
  
}

#' Przedział ufności
#' Funkcja \code{Confidence.interval95} generuje przedział ufności 
#' 
#' @param price.vector wektor cen opcji
#' 
#' @return Wartość minimalna i maksymalna przedziału ufności
#' @export
Confidence.interval95 <- function(price.vector){
  
  price.sd  <- sd(price.vector)
  price.mean <- mean(price.vector)
  M    <- length(price.vector)
  min.interval <- price.mean - 1.96 * price.sd /sqrt(M)
  max.interval <- price.mean + 1.96 * price.sd /sqrt(M)
  interval     <- c(min.interval, max.interval)
  names(interval) <- c("Min.confidence.interval","Max.confidence.interval")
  interval
}

#' Wycena opcji Call w Modelu Hestona
#' Funkcja \code{Heston.Monte.Carlo.call.price.option} wyznacza cenę opcji call w modelu Hestona
#' 
#' @param M liczba 
#' @param n liczba trajektorii procesu Wienera
#' @param S0 S0
#' @param V0 V0
#' @param r  stopa procentowa
#' @param sigma zmienność
#' @param a  param b modelu
#' @param b  param a modelu
#' @param rho korelacja procesów Wieniera
#' @param Time czas
#' @param K cena wykonania
#' @export
Heston.Monte.Carlo.call.price.option <- function(M, n, S0 = 100, V0 = 0.09, r = 0.05, sigma = 1.0, a = 2, b=0.09, rho = -0.3, Time = 5.0, K = 100){
  
  stopifnot(M>=10^3,n>=400, rho <= 1, rho >= -1)
  
  
  delta.t <- Time / n
  discount.factor <- exp( - r * Time)
  alpha1 <- a * delta.t
  alpha2 <- sigma * sqrt(delta.t)
  theta1 <- 1 + r * delta.t
  
  x<- rep(n,M)
  Q1<-lapply(x, Heston.model.one.iteration, S0 = S0, V0 = V0,  sigma = sigma, rho = rho, K = K, b=b, delta.t = delta.t, alpha1 = alpha1 ,alpha2= alpha2, theta1= theta1 ,df = discount.factor)
  Q <- unlist(Q1)
  
  value <-mean(Q)
  names(value) <- "Heston.call.price"
  interval <- Confidence.interval95(Q)
  c(value,interval)
  
}


#' Symulacje Monte Carlo - model Hestona 
#' Funkcja \code{Heston_call_MC} oblicza cenę z europejskiej opcji call w modelu Hestona
#' @param file.name nazwa pliku tekstowego z danymi
#' 
#' @return cena opcji call w modelu hestona
#' @export

Heston_call_MC <- function(file.name){
  
  file <- paste0(file.name,".txt")
  data <-read.csv(file, sep = ";")
  stopifnot(nrow(data) ==1)
  stopifnot(c("M","n", "S0","V0", "r","T","K","sigma","a","b","rho") %in% colnames(data))
  M     <- data$M
  n     <- data$n
  S0    <- data$S0
  V0    <- data$V0
  r     <- data$r
  sigma <- data$sigma
  a     <- data$a
  b     <- data$b
  rho   <- data$rho
  Time  <- data$T
  K     <- data$K
  
  heston   <- Heston.Monte.Carlo.call.price.option(M=M, n = n, S0 =S0, V0 =V0, r = 0.05, sigma = sigma, a = a , b = b , rho = rho, Time = Time, K = K)
  names(M) <- "Number.of.symulations"
  names(n) <- "Trajectory.of.Brownian.motion"
  c(M,n,heston)
  
}
############################################################
######## Test funkcji Heston_call_MC #######################
############################################################
Heston_call_MC("CW3_data")