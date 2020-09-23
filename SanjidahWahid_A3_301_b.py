'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 3 - #301 part b

Integrals in Our Solar System: 
Compute the distance Halley’s comet travels in one complete orbit around the Sun
'''

import numpy as np
import matplotlib.pyplot as plt

a_2 = 2693165649421 # From 101; in meters
b_2 = 682178399744 # From 101; in meters

M = int(input("Number of partitions for trapezoidal evaluation M is: "))

def function(theta):
    return np.sqrt((a_2 * np.sin(theta))**2 + (b_2 * np.cos(theta))**2)

#test function
# def function(t):
#     return 3*t**2 + 1

def dfunction(t):
    return 6*t

# def dfunction(theta):
#     return -((b_2-a_2)*(b_2+a_2)*np.cos(theta)*np.sin(theta))/np.sqrt(a_2**2*np.sin(theta)**2 + b_2**2*np.cos(theta)**2)

test = np.linspace(0,np.pi*2, 100)

plt.plot(test,function(test))

a = 0. #alpha
b = 2*np.pi #beta
g = (b - a)/M

#Trapezoidal rule evalation
# def trap(M):
#     I_trap_end = 1/2*(function(0)+function(2*np.pi))
#     I_trap_val = 0.
#     for n in range(1,M):
#         I_trap_val += function(alpha+n*g)
#     I_trap = g * (I_trap_val + I_trap_end)
#     return I_trap

                    
def trap(M):
    end = 1/2 * (function(a) + function(b))
    odd = 0
    even = 0
    d = (b - a)/M
    for n in range(1,M//2):
        odd += function(a + (2*n-1)*d)
        even += function(a + 2*n*d)
    odd += function(b - d)
    return g* (end + odd + even)
    
exact = 8*(np.pi)**3 + 2*np.pi
print(str(exact - trap(M)))

print(str(trap(M)) + " meters in length with trapezoidal")

def mod_I_trap(M):
    return trap(M) + g**2/12 * (dfunction(a_2)-dfunction(b_2))

def abs_err(M):
    return 1/3 * (np.abs(trap(M)-trap(2*M)))

print(abs_err(int(2*M)))
print(abs_err(M))

delta = 1
epsilon = 1e-5

while delta >= epsilon:
 #   q = trap_eval(M)
#    relative = np.abs(np.abs(trap_eval(int(M/2))-q)/q)
    relative = np.abs(abs_err(M))/np.abs(trap(2*M))
    M *= 2
    print(str(M) + " :M")
    print(str(relative) + " :delta")
    delta = relative

print(str(delta) + " :final delta")

# C = 1e-16
# manual_abs_err = C * ((b - a)/(M*2))**2
# print(abs_err_calc(M))
# print(manual_abs_err)

# err_est = 1/3 * (mod_I_trap(M) - mod_I_trap(2*M))
# print(str(err_est) + " absolute error using mod trapezoidal rule")

# print(str(mod_I_trap(M)) + " meters in using mod trapezoidal rule")

# err_est_mod = 1/3* (mod_I_trap(M) - mod_I_trap(2*M))
# print(str(err_est_mod))

#rel_err = abs(exact_int - mod_I_trap(M)) / exact_int
#print(rel_err)
