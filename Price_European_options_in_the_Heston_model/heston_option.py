# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#n = 1000
rho = -0.3
Time = 5.0
S0 = 100
V0 = 0.09
r = 0.05
sigma = 1.0
a = 2
b = 0.09
rho = -0.3
Time = 5.0
K = 100
import numpy as np
#import  scipy.stats as stats
def Heston_one_iteration(n , S0 = 100 , V0 = 0.09, r = 0.05 , sigma = 1.0 , a = 2, b = 0.09, rho = -0.3  , Time = 5.0, K = 100):

    #  Generowanie zmiennych normalnych
    Zv = np.random.normal(size = n)
    Z2 = np.random.normal(size = n)
    Zs = rho * Zv + np.sqrt(1 - rho ** 2) * Z2
    
    # Trajektorie Vt - schemat Mildsteina
    deltat = Time / n
    V = np.zeros(n)
    V[0] = V0 
    beta = (sigma ** 2 )* deltat * (np.power(Zv,2) -1) * 0.25
    alpha1 = a * deltat
    alpha2 = sigma * np.sqrt(deltat)
    theta1 = 1 + r * deltat
    Zsprime = Zs * np.sqrt(deltat)
    for i in range(1, n):
        V[i] = V[i-1] + alpha1 * (b - max(V[i-1],0))+ alpha2 * np.sqrt(max(V[i-1],0)) *Zv[i -1] + beta[i-1]

    # Wyznaczenie  wartosci procesu S w chwili T - schemat Eulera
    Vprime2 = np.where(V < 0 ,0,V)
    S = S0 * np.prod(theta1 + np.sqrt(Vprime2) * Zsprime)
    ## Wyznaczanie wartosci opcji kupna 
    C = np.exp( - r * Time) * max(S - K , 0)
    return C

def Confidence_interval95(price_vector):
    price_sd = np.std(price_vector)
    price_mean = np.mean(price_vector)
    M = len(price_vector)
    max_interval = price_mean + 1.95 * price_sd  / np.sqrt(M)
    min_interval = price_mean - 1.95 * price_sd / np.sqrt(M)
    interval = {'min_confidence_interval': min_interval, 'max_confidence_interval': max_interval}
    return  interval

def Heston_Monte_Carlo_call(M, n , S0 = 100 , V0 = 0.09, r = 0.05 , sigma = 1.0 , a = 2, b = 0.09, rho = -0.3  , Time = 5.0, K = 100):
    prime = np.zeros(M)
    for i in range(M):
        prime[i] = Heston_one_iteration(n)    
    interval = Confidence_interval95(prime)
    interval["heston_call_price"] = np.mean(prime)
    for s in interval:
        print(s,":",interval[s])   
 
my_data =  np.genfromtxt("CW3_data.txt", delimiter = ';', names = True)

M = my_data['M']
n = my_data['n']
a = my_data['a']   