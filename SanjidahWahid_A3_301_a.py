'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 3 - #301 part a

Integrals in Our Solar System
'''

import numpy as np

sun_radius = 695500
earth_radius = 6371
earth_orbit_rad = 149.6e6
G = 6.674e-11
earth_mass = 5.9722e24
sun_mass = 1.9884e30

freefall_distance = earth_orbit_rad - (sun_radius + earth_radius)

N = 80
a = earth_radius
b = freefall_distance + earth_radius

E_m = G * earth_mass * sun_mass * (1/earth_radius - 1/(freefall_distance + earth_radius))

def internal(x):
    return ((2*x) / (G*sun_mass)) *  np.sqrt((2/earth_mass) * (E_m + (G*sun_mass*earth_mass)/x))

def external(x):
    return ((-x**2) / (G*sun_mass)) *  np.sqrt((2/earth_mass) * (E_m + (G*sun_mass*earth_mass)/x))

h = (b-a) / N

I_S_end = 1/3 * (internal(a) + internal(b))

def I_Simp():
    Simp = 0
    for i in range(1,N):
        if i%2 == 0:
            Simp += 2/3 * internal(a + i*h)
        else:
            Simp += 4/3 + internal(a + i*h)
    return Simp

I_S = h * (I_Simp() + I_S_end)

delta_t = (external(b) - external(a)) + I_S

print("delta_t = " + str(delta_t))









# the final answer should be based on the equation you got from the integration 
# by parts, it would be any constant out front times the sum of the 'uv' term 
# plus the integral in your f(x)
#### print(leading constant * ( uv term + integral))

# def Isimp(f, a, b, N=10000):
#     h = (b-a)/N
#     S = f(a) + f(b) + 4*f(a + h*(N-1))
#     for k in range(1,N//2):
#         S += 4*f(a + h*(2*k-1)) + 2*f(a + h*2*k)
#     return h/3*S
# print(Isimp)
