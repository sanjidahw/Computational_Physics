'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 3 - #301 part a

Integrals in Our Solar System
'''

from numpy import sqrt

sun_radius = 695500e3 # meters
earth_radius = 6371e3 # meters
earth_orbit_rad = 149.6e9 
G = 6.674e-11
M = 1.9891e30
earth_mass = 5.9722e24 # kg
sun_mass = 1.9884e30 # kg

N = 80
xf = sun_radius + earth_radius
xi = earth_orbit_rad

#Let's solve for the integral first
def f(x):
    return sqrt((xi-x)/x)

# Simpson's rule:
    
h = (xf-xi)/N

s = (1/3) * (f(xi) + f(xf))  #This is for the end points, we treat them different

s_odd = 0
s_even = 0

for n in range(1,N//2):
    s_odd += f(xi + (2*n-1)*h)
    s_even += f(xi + 2*n*h)
    
s_odd += f(xf-h)

s_answer = h*(s + (4/3)*s_odd + (2/3)*s_even) 
#print(s_answer) #value of the integral

# Now let's find the solution for the full expression
I_S = (-sqrt(xi/(2*G*M)) * (-2*sqrt(xf*(xi - xf)) + s_answer)) 

days = I_S / (60*60*24)
print('Simpson with N = {:4} steps, I = {:.0f}, days = {:.0f} days'.format(N,I_S,days))
                           
