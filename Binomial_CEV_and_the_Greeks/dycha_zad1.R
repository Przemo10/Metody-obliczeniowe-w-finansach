#' Process X
#' Funkcja \code{xtransform.process} wyznacza drzewo początkowe procesu X
#' 
#' @param M liczba podziałów 
#' @param S0
#' @param beta
#' @param sigma
#' @param Time
#' 
#' @return Proces X  z algorytmu modelu CEV
#' 
xtransform.process <- function(M,S0 = 105, beta =1, sigma =0.25, Time = 1) {
  
  alpha  <- 1 - beta/2
  delta   = sigma * S0^(alpha)
  X0     <- S0^alpha/(alpha * delta)
  delta.t = Time/M
  shift   = sqrt(delta.t)
  shift2  = 2 * shift
  xtreeprocess <- c()
  #xtreeprocess[1] <- X0
  i <- 1
  while( i <=(M+3)){
    
    len <- length(xtreeprocess)
    shiftvalue <- (i -1) * shift
    xtreeprocess[(len+1) :(len+i)]  <- seq.int(X0 - shiftvalue, X0 + shiftvalue,shift2)
    i <- i + 1
  }
  
  xtreeprocess
  
}
#' Proces odwrotny do procesu X
#' Funkcja \code{inverse.transform} wyznacza "proces odwrotny" do procesu X
#' 

inverse.transform <- function(transform.process, alpha = 1/2, delta =  2.561738){
  
  # stopifnot(is.vector(transform.process) && is.numeric(transform.process))
  stopifnot(alpha < 1 && alpha > 0)
  
  inverse.transform <- ifelse(transform.process > 0, (transform.process * alpha * delta)^(1/alpha), 0)
  
  
  
  inverse.transform
}

#' Suma od 1 do n
#' Funkcja \code{sum.1tok.plus} oblicza sumę liczb od 1 do n + k
#' @param k liczba do której należy sumować
#' @param l parametr przesunięcia
#' 
#' @return sumę od 1 do n + k
#' 
sum.1tok.plusl <- function(k,l =0){
  
  k*(k+1)/2 + l
  
}

#' Funkcja wyznaczająca miarę martyngałową w modelu CEV
#' 
#' @param sprocess drzewo procesu S
#' @param N liczba podziałów drzewa
#' @param r stopa procentowa modelu B-S
#' @param Time czas
#' 
#' Funkcja wykorzystuje charakterystyczną własność drzewa w modelu CEV.
#' Prawdopodobieństwa martyngałowe są takie same w zależności od parzystości paramteru N
#' W związku z tym wystarczy wyznaczyć prawdopodobieństwa na dwóch ostatnich okresach.
#' 
#' @return Wektor prawdopodobieństw rozmiaru N-1
#' 

probability.risk <-  function(sprocess, N, r = 0.09, Time = 1){
  
  deltat =  Time / N
  
  stopifnot( sum.1tok.plusl(N+3) == length(sprocess) || sum.1tok.plusl(N+4) == length(sprocess))
  
  
  s <-  sprocess[sum.1tok.plusl(N+2,1):sum.1tok.plusl(N +3)]
  s1 <- sprocess[sum.1tok.plusl(N+1,1):sum.1tok.plusl(N +2)]
  
  diff = (s[-1] -s[-(N+3)]) 
  
  
  prob<-  ifelse(diff != 0,  (s1 * exp(r * deltat) - s[-(N + 3)]) /
                   (diff) , 0)
  
  prob <- ifelse(prob <  1, prob, 1)
  prob <- ifelse(prob >= 0, prob, 0)
  
  prob
}

#' 
#' Funkcja \code{calc.option.value} wyznacza  wartości opcji dla j-tego okresu
#' @param j j-ty okres drzewa
#' @param vec wektor cen opcji
#' @param prob wektor prawdopodobieństw
#' @param DF czynnik dyskątowy
#' 
#' @return Wartości opcji dla poziomu j
#' 
calc.option.value <- function(j,vec, prob , DF){
  
  
  value <- (vec[-1] * prob + vec[-(j+1)] * (1-prob) ) * DF
  value
  
}
#' Funckja \code{cut.prob.vector} skraca wektor prawdopodobieństw w celu zachowania jednakowych długości.
#' @param prob  wektor prawdopodobieństw
#' @param j i-ta pozycja wektora
#' 
#' Funkcja skraca wektor prawdopodobieństw o j-tą oraz pierwszą pozycję wektora.
#' 
#' @return Skrócony wektor prob
cut.prob.vector <- function(prob, j){
  
  prob <- prob[-j]
  prob <- prob[-1]
  prob
  
}


#' Grecki parametr Delta
#' 
#' Funkcja \code{Greeks.Delta} wyznacza grecki parametr Delta.
#' 
#' @param sprocess wektor procesu S
#' @param option.value wector wartości opcji V
#'  
#' @return grecki parametr Delta
#
Greeks.Delta <- function(sprocess,option.value){
  
  delta <- (option.value[6] - option.value[4])/(sprocess[6]- sprocess[4])
  
  delta  
  
  
}

#' Grecki parametr Theta
#' 
#' Funkcja \code{Greeks.Theta} wyznacza grecki parametr Theta.
#' 
#' @param option.value wektor wartości opcji
#' @param deltat przedział czasowy
#'   
#' @return grecki parametr Theta
#' 
Greeks.Theta <- function(option.value, deltat){
  
  theta <- (option.value[13] - option.value[1])/( 4 * deltat)
  theta 
}


#' Grecki parametr Gamma
#' 
#' Funkcja \code{Greeks.Gamma} wyznacza grecki parametr Gamma.
#' 
#' @param sprocess wektor procesu S
#' @param option.value wector wartości opcji V
#'  
#' @return grecki parametr Gamma
#
Greeks.Gamma <- function(sprocess, option.value){
  
  v1 <- (option.value[6] - option.value[5])/(sprocess[6]- sprocess[5])
  v2 <- (option.value[5] - option.value[4])/(sprocess[5]- sprocess[4])
  t  <- (sprocess[6]- sprocess[4])
  gamma <- 2*(v1 - v2)/t
  gamma 
}


#' Greckie parametry oraz wartość opcji w modelu CEV
#' 
#' Funkcja \code{call.cev.option} wyznacza Greckie parametry oraz wartość opcji 
#' 
#' @param M liczba okresów
#' @param S0 
#' @param r
#' @param sigma
#' @param Time
#' @param K
#' @param beta
#' @param print.option.value zaznaczyć true jeśli chcemy wyznaczyć zobaczyć całę drzewo opjci
#' @param print.xprocess czy wyświetlić proces X
#' @param print.sprocess czy wyświetlić proces S
#' @return M,V0, Delta, Theta, Gamma
#' 
call.cev.option <-
  function(M,
           S0 = 105,
           r = 0.09,
           sigma = 0.25 ,
           Time = 1 ,
           K = 100,
           beta = 1,
           print.option.value = FALSE,
           print.xprocess = FALSE,
           print.sprocess = FALSE) {
    
  stopifnot(beta <= 2 , beta >= 0)
  stopifnot(is.numeric(S0),is.numeric(sigma),is.numeric(Time),Time >= 0, K >0)
  
  alpha1  <- 1 - beta/2
  delta1 = sigma * S0^(alpha1)
  X0 <- S0^alpha1/(alpha1 * delta1)
  delta.t <- Time / M
  
  discount.factor <- exp(- r * delta.t)
  xprocess <- xtransform.process(M,S0, beta, sigma, Time)
  sprocess <- inverse.transform(xprocess, alpha =  alpha1,delta = delta1)
  
  p  <-  probability.risk(sprocess, N = M  , r = r, Time = Time)
  p1 <-  probability.risk(sprocess, N = M-1, r = r, Time = Time)
  
  
  option.value <- vector( , sum.1tok.plusl(M+3))
  S  <- sprocess[sum.1tok.plusl(M+2,1):sum.1tok.plusl(M +3)]
  option.value[sum.1tok.plusl(M+2,1):sum.1tok.plusl(M +3)] <- ifelse(S - K > 0, S-K,0)
  
  i <- M+2
  while (i >1){
    
    v<-option.value[sum.1tok.plusl(i,1):sum.1tok.plusl(i+1)]
    
    if( ((M-i) %% 2) == 0){
      value <- calc.option.value(j =i,vec = v,prob = p,DF = discount.factor)
      p <- cut.prob.vector(prob = p,j = i)
    } else {
      value <- calc.option.value(j =i,vec = v,prob = p1,DF = discount.factor)
      p1 <- cut.prob.vector(prob = p1, j =i)
    }
    option.value[sum.1tok.plusl(i-1,1):sum.1tok.plusl(i)] <- value
    
    i <- i-1
  }
  
  if( ((M-1) %% 2) == 0 ){
    value <- (p *  option.value[3] + ( 1- p) * option.value[2]) * discount.factor
    
  } else{
    value <- (p1 *  option.value[3] + ( 1- p1) * option.value[2]) * discount.factor
  }
  
  option.value[1] <- value
  

  Delta <- Greeks.Delta(sprocess, option.value)
  Gamma <- Greeks.Gamma(sprocess, option.value)
  Theta <- Greeks.Theta(option.value, deltat = delta.t)
 
  ### Wyswietlanie   
  if(print.xprocess){
    
    print(paste0("Drzewo procesu X"))
    print(xprocess)
  }
  
  if(print.sprocess){
    print(paste0("Drzewo procesu S"))
    print(sprocess)
  }
  if( print.option.value ){
    print(paste0("Wartosci opcji"))
    print(option.value)
  }
  Delta <- round(Delta,6)
  Gamma <- round(Gamma,6)
  Theta<- round(Theta,6)
  result.vector<-c(M,option.value[5],Delta, Theta, Gamma)  
  names(result.vector) <- c("M","V0", "Delta", "Theta", "Gamma")
  result.vector
}

#' Funkcja \code{binomial_CEV} wczytuje dane z pliku tekstowego i wywołuje funkcję \code{call.call.option}
#' 
#' @param filename nazwa pliku (bez rozszerzenia)
#' 
#' Plik tesktowy oddzielany średnikiem
#' Wymagane kolumny:
#' M, S0, r, T, K, beta
#' 
#' @return liczba M, V, Delta, Theta , Gamma
#' 
#' 
binomial_CEV <-
  function(file.name){
    
    file <- paste0(file.name,".txt")
    data <-read.csv(file, sep = ";")
    stopifnot(nrow(data) ==1)
    stopifnot(c("M","S0","r","T","K","beta") %in% colnames(data))
    N     <- data$M
    S0    <- data$S0
    r     <- data$r
    sigma <- data$sigma
    Time  <- data$T
    K     <- data$K
    beta  <- data$beta
    
   call.cev.option( M     = N,
                    S0    = S0,
                    r     = r,
                    sigma = sigma,
                    Time  = Time,
                    K     = K,
                    beta  = beta)
  
  
  }

#'@example
binomial_CEV("CW1_data")