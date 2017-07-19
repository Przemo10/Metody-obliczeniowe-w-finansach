#' Wartość europejskiej opcji Put metodą Monte Carlo.
#' Funkcja \code{PutMC} wyznacza wartość opcji put dla metod Monte Carlo.
#'  
#' @param N liczba symulacji
#' @param S0 wartość opcji w chwili 0
#' @param sigma zmienność
#' @param r stopa procentowa
#' @param Time czas
#' @param QMC parametr dodatkowy, czy uwzględnić przy symulacji metodę QMC  
#' @param confidence.interval.MC parametr dodatkowy, czy zwrócić informację o przedziale ufności 95% 
#' 
#' @return wartość opcji put oraz błąd  oraz dodatkowo przedzial ufnosci metody MC
#' 
#' @examples 
#' PutMC(100, QMC = FALSE)
#' PutMC(100, confidence.interval.MC = TRUE ) 
#' @export
PutMC<- function(N ,S0 = 50, K = 50,sigma= 0.3,r =0.05,Time = 0.5, QMC = TRUE, confidence.interval.MC = FALSE){
  
  stopifnot(QMC == TRUE | QMC == FALSE)
  stopifnot(confidence.interval.MC == TRUE || confidence.interval.MC == FALSE)
  stopifnot(N > 1)
  
  a  <-  r - sigma ^ 2 / 2
  b  <-  sigma * sqrt(Time)
  df <-  exp(-r * Time)
  
  bs.price <- Black.Scholes.put.formula(S0,K,r,Time,sigma)
  
  u1  <- runif(floor(N))
  len <- (length(u1))
  U   <- Box.Muller(u1[1:(len / 2)], u1[(len / 2 + 1):len])
  X   <- Put.price.vector(U, S0, K, Time, a, b, df)
  
  put.price <- mean(X)
  names(put.price)[1] <- "Put.price.runif"
  error.runif <-  Error.monte.carlo(put.price, bs.price)
  names(error.runif) <- "error.runif"
  put.price <- append(put.price,error.runif)
  
  if(confidence.interval.MC){
    conf <- Confidence.interval95(X)
    put.price <- append(put.price, conf)
  }
  
  if(QMC){
    h1 <- Halton.sequence(floor(N / 2), prime = 2)
    h2 <- Halton.sequence(floor(N / 2), prime = 3)
    H  <- Box.Muller(h1, h2)
    X  <- Put.price.vector(H, S0, K, Time, a, b, df)
    halton.put.price <- mean(X)
    names(halton.put.price) <- "Put.price.halton"
    error.halton <- Error.monte.carlo(halton.put.price, bs.price)
    names(error.halton) <- "error.halton"
    halton.put.price <- append(halton.put.price, error.halton)
    put.price <- append(halton.put.price, put.price)
  }
  
  put.price
}

#' Wycena Opcji Put w Modelu Blacka Scholesa
#' Funkcja \code{Black.Scholes.put.formula} wyznacza wartość opcji w modelu B-S
#' 
#' @param S0
#' @param K
#' @param r
#' @param Time
#' @param sigma 
#' 
#' @return C0  
#' @export    
Black.Scholes.put.formula <- function(S0 = 50, K = 50 ,r = 0.05, Time  = 0.5 , sigma = 0.3 ){
  
  d1 <- (log(S0/K) + (r  + sigma^2/2)*Time)/(sigma * sqrt(Time))
  d2 <-  d1 - sigma * sqrt(Time) 
  Put <- K * exp(-r * Time) * pnorm(-d2) - S0 * pnorm(-d1)
  names(Put) <- "BS.put"
  Put
}

#' Generowanie zmiennych N(0,1) z rozkładu jednostajnego
#' 
#' Funkcja \code{Box.Muller} generuje zmienne z rozkładu normalnego
#' 
#' @param unif1 wektor dł k z rozkładu jednostajnego
#' @param unif2 wektor dł k z rozkładu jednostajnego
#' 
#' @return wektor długości 2k z rozkładu normalnego
#' 
#' @export
Box.Muller <- function(unif1, unif2){
  
  stopifnot(length(unif1) == length(unif2))
  
  Theta <-  2 * pi * unif1
  R     <-  sqrt(- 2 * log(unif2) )
  X     <-  R * cos(Theta)
  Y     <-  R * sin(Theta)
  c(X,Y)
}

#' Generowanie wektora cen opcji put dla zadanego wektora z rozkładu normalnego
#' \code{Put.price.vector}
#' @param norm.dist.vector wektor rozkładu normalnego
#' @param S0 
#' @param K cena wykonania
#' @param Time
#' @param a = (r - sigma ^ 2 / 2)
#' @param b =  sigma * sqrt(Time)
#' @param df czynnik dyskontowy (exp(-r * Time))
#' 
#' @return Put.price.vector
#' @export
Put.price.vector <- function(normal.dist.vector,S0 = 50,K = 50,Time = 0.5 ,a = 0.05,b = 0.212132, df = 0.9753099){
  
  St <- K  - S0 * exp(a * Time + b * normal.dist.vector)
  St <- ifelse(St >0, St,0)
  P0 <- df * St
  P0
}

#' Genorowanie liczby z ciągu Haltona
#' Funkcja \code{Halton.number} generuje n-ty wyraz ciągu Haltona dla ustalonej liczby pierwszej p
#' Sposób implementacji z \link{https://en.wikipedia.org/wiki/Halton_sequence}
#' @param n 
#' @param p liczba pierwsza
#' 
#' @return  liczba Haltona 
#' @export
Halton.number <- function(n,p){
  
  f <- 1
  r <- 0
  while(n > 0){
    f <- f/p
    r <- r + f * (n %% p)
    n <- floor(n/p)
    
  }
  
  r
}

#' Ciąg Haltona
#' Funkcja \code{Halton.sequence} generuje ciąg haltona długości k  dla ustalonej liczby pierwszej
#' @param k długość ciągu
#' @param prime liczba pierwsza
#' 
#' @return ciąg Haltona
#' 
#' @export
Halton.sequence <- function(k,prime){
 
  sapply(1:k, Halton.number, p = prime)
  
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


#' Błąd wyceny Monte Carlo
#' 
#' Funkcja \code{Error.monte.carlo} zwraca błąd metody monte carlo
#' 
#' @param monte.carlo.estimate cena Monte Carlo
#' @param black.sholes.price   wycena formuł B-S
#' 
#' @return abs(monte.carlo.estimate - black.sholes.price)/black.sholes.price
#' @export
Error.monte.carlo <- function(monte.carlo.estimate, black.sholes.price){
  
   abs(monte.carlo.estimate - black.sholes.price)/black.sholes.price
  
}
################PROGRAM######################################
FILE.NAME <- "CW2_data.txt"
data <-read.csv(FILE.NAME, sep = ";")
stopifnot(nrow(data) ==1)
stopifnot(c("M","S0","r","T","K") %in% colnames(data))
N     <- data$M
S0    <- data$S0
r     <- data$r
sigma <- data$sigma
Time  <- data$T
K     <- data$K
if ("QMC" %in% colnames(data)) {
  qmc <- data$QMC
} else{
  qmc <- TRUE
}
if ("confidence.interval.MC" %in% colnames(data)) {
  conf <- data$confidence.interval.MC
} else{
  conf <- FALSE
}


put.mc   <- PutMC(N, S0 = S0, K = K, sigma = sigma, r = r, Time = Time, QMC = qmc, confidence.interval.MC = conf)
names(N) <- "number.of.symulations"
black.scholes.put <- Black.Scholes.put.formula(S0,K,r,Time,sigma)
result <- c(N,black.scholes.put, put.mc)
options(scipen = 999)
result
