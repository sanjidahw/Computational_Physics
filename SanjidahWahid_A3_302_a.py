'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 3 - #302 part a and b

Madelung Constant from Integratio
'''

import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt

#part (a)

mp.dps = 16 
mp.pretty = True  

Ngrid = 1000

z = np.linspace(0.000001,0.999999,Ngrid)

theta4 = np.vectorize(lambda qval: mp.jtheta(4, 0, qval))   # n = 4, z = 0

dztheta4 = np.vectorize(lambda qval: mp.jtheta(4, 0, qval, 1))

qdqtheta4 = np.vectorize(lambda qval: -0.25*mp.jtheta(4,0,qval,2))

def f(z):
    s = np.sqrt(z) / (1-np.sqrt(z))
    return 3*theta4(np.exp(-s))**2 * qdqtheta4((np.exp(-s))) * np.sqrt(s/np.pi) \
           / (np.sqrt(z)*(1 - np.sqrt(z))**2)

plt.plot(z,f(z),"-r")
plt.show()


##############################################################################
##############################################################################


#part (b) 

# Gaussian quadrature from 0 to 1
def gaussxw(N):
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))
    
    #Find roots with Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)* (p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))
        
    #  weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
    
    return x,w

# Gaussian quadrature from a to b
def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a), 0.5*(b-a)*w

# Madelung constant value for reference
M = -1.74756459463318
a = 0.
b = 1.
N = int(input("Select an integer value for N: "))

x,w = gaussxwab(20,a,b)
I = np.sum(w*f(x))

I_G = np.sum(w*f(x))
print(str(I_G) + " : Gaussian quadrature evaluation for N = " + str(N))

plt.plot(x,w*10,"xb")

M_bar_30 = -1.7286350779651571
M_bar_100 = -1.7418198158362812

err_g = np.abs(np.abs(M - I_G)/M) * 100
print(str(err_g) + " error with Gauss's rule for N = " + str(20))
err_bar_30 = np.abs(np.abs(M_bar_30 - I_G)/M_bar_30) * 100
print(str(err_bar_30) + " error with partial sums evaluation for N = 30")
err_bar_100 = np.abs(np.abs(M_bar_100 - I_G)/M_bar_100) * 100
print(str(err_bar_100) + " error with partial sums evaluation for N = 100")

a = 0.0000001
b = 0.9999999
h = (b-a)/N

I_S_end = 1/3 * (f(a) + f(b)) #endpoints

def I_S():
    I_S_val = 0.
    for i in range(1,N):
        if i%2 == 0:
            I_S_val += 2/3 * f(a + i*h)
        else:
            I_S_val += 4/3 * f(a + i*h)
    return I_S_val

I_simp = h*(I_S() + I_S_end)

err = np.abs(np.abs(M - I_simp) / M) *100
print(str(err) + " error with Simpson's rule for N = " + str(N))