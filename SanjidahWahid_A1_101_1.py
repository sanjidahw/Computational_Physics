'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 1 - #101

Newman, Exercise 2.6: Planetary Orbits with the following modifications:
    (a). While you must perform the calculation of v2, do not turn it in.  Instead, demonstrate thatyour solution is valid by verifying that it conserves mechanical energy, see below.
    (b). Add  two  items  to  the  output:  mechanical  energy  per  mass E/m,  at  perihelion  and  at aphelion.  Equality of the two verifies the solution from (a).
    (c). Make sure the period is output in units of years
'''

from numpy import sqrt, pi

G = 6.6738*(10**-11)    #Newton's gravitational constamt
M = 1.9819*(10**30)    #Mass of the Sun

l1 = float(input("Enter the distance to the Sun in meters (This is l1): "))
v1 = float(input("Enter the velocity at perihilion in meters per second (This is v1): "))

v2 = ((2*G*M)/(v1*l1)-sqrt(((2*G*M)/(v1*l1))**2-4*((2*G*M)/l1-v1**2)))/2
l2 = l1*v1/v2
a = 1/2*(l1+l2)     #Semi-major axis
b = sqrt(l1*l2)     #Semi-minor axis
T = (2*pi*a*b)/(l1*v1)*((1/365)*(1/24)*(1/60)*(1/60))    #Orbital period converted to years
e = (l2-l1)/(l2+l1)     #Orbital eccentricity
E1 = (1/2)*v1**2 - G*M/l1    #Total Energy with v1
E2 = (1/2)*v2**2 - G*M/l2    #Total Energy with v2

print('Part B: Mechanical energy per unit mass at perihilion = {:.4E} J'.format(E1))
print('Part B: Mechanical energy per unit mass at apihilion = {:.4E} J'.format(E2))

print(f'Part C: l2 = {l2} m\n v2 = {v2} m/s\n T = {T} years\n e = {e}\n')

