'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 5 - #501 

Newman, Exercise 6.16: The Lagrange Point
'''

import numpy as np
import matplotlib.pyplot as plt

G = 6.674e-11 
M = 5.974e24 
m = 7.348e22 
R = 3.844e8 
w = 2.662e-6


#Newton's Method

print("Newton's Method in meters:")
def f(r):
	return -w**2*r**5 + 2*w**2*R*r**4 - w**2*R**2*r**3 + G*r**2*(M -m) - 2*r*R*G*M +G*M*R**2
	
def df(r):
	return -5*w**2*r**4 + 8*w**2*R*r**3 - 3*w**2*R**2*r**2 + 2*G*r*(M - m) - 2*R*G*M
	
del_N = 1e-8
epsilon_N = 1.0e-4
rguess = 3.2e8 # meters

while epsilon_N>del_N:
	fr = f(rguess)
	dfr = df(rguess)
	dr = - fr/dfr
	rguess += dr
	epsilon_N = np.abs(dr)
	print(rguess)
	

# Part a plotting
fig_a = plt.figure(figsize = (12,12))
r = np.linspace(0.1596003 * R, R, 100) # L1 distance
f_vals = list(map(f, r))
plt.grid()
plt.plot(r, f_vals, "--", color = "lightseagreen")
idx = np.argwhere(np.diff(np.sign(f_vals))).flatten()
plt.plot(r[idx], 0, 'o', color = "black")
plt.title("Part a: The Lagrange point of L1 satisfies f(r)", fontsize = 20)
plt.xlabel("Distance from the center of the Earth (m)", fontsize = 18)
plt.ylabel("f(r)", fontsize = 18)
plt.show()
print(type(idx))	



#Relaxation Method

print("Relaxation Method in meters:")

def F(r):
	#return (G*M/r**2) - G*m/(R-r)**2 - r*w**2
	return (w**2*r**5 - 2*w**2*R*r**4 + w**2*R**2*r**3 + G*r**2*(m - M) - G*M*R**2)/(-2*R*G*M)

delta = 1e-8
epsilon = 1.0

r = 3.2e8

while epsilon>delta:
	rnew = F(r)
	epsilon = np.abs(r/rnew - 1)
	r = rnew
	print(np.abs(r))



#Bisection Method

print("Bisection Method in meters:")
def func(r):

	return (G*M/r**2) - G*m/(R-r)**2 - r*w**2
	
a = -3.0e5
b = 3.2e10

def findroot(func, a, b, del_B = 1e-8):
	 epsilon = 1.0
	 while epsilon>del_B:
		 mid = 0.5*(a + b)
		 if func(mid)/ func(a)<0:
			 b = mid
		 else:
			 a = mid
		 epsilon = 2*np.abs(a - b)/(np.abs(a) + np.abs(b))
	 return 0.5*(a + b)

print(findroot(func, a, b))
