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

G = 6.6738*(10**-11)    # Newton's gravitational constamt
M = 1.9819*(10**30)    # Mass of the Sun
## can give in scientific notation, e.g. 1.9891e30, else python has to compute
## actually is 91 rather than 19


l1 = float(input("Enter the distance to the Sun in meters (This is l1): "))
v1 = float(input("Enter the velocity at perihilion in meters per second (This is v1): "))

# Extract v2 using quadratic formula
a=1
b=-((2*G*M)/(v1*l1))
c=-(v1**2-((2*G*M)/l1))
## superfluous parenthesis make code harder to read


# Find the roots
D=b*b-4*a*c
root1=(-b-D**.5)/2
root2=(-b+D**.5)/2
v2=min([root1,root2])

l2=(v1*l1)/v2    # Farthest point of orbit

E1 = (1/2)*v1**2 - G*M/l1    # Total Energy with v1 (perihilion)
E2 = (1/2)*v2**2 - G*M/l2    # Total Energy with v2 (apihilion)

a=1.0/2*(l1+l2)    # Semi-major axis
b=sqrt(l1*l2)    # Semi-minor axis
T=(2*pi*a*b)/(l1*v1)/(24*60*60)/365    #Orbital period converted to years
## should be sidereal year rather than calendar year

e=(l2-l1)/(l2+l1) #Orbital eccentricity

print("------------Part (b) and Part (c)------------")

print('Part B: Mechanical energy per unit mass E/m at perihilion = {:.4E} J'.format(E1))
print('Part B: Mechanical energy per unit mass E/m at apihilion = {:.4E} J'.format(E2))
print('Part C: Velocity at aphelion (v2) = {:.4e} m/s'.format(v2))
print('Part C: Farthest point of orbit (l2) = {:.4e} m'.format(l2))
print('Part C: Time period (T) = {:.4e} years'.format(T))
print('Part C: Orbital eccentricity (e) = {:.6e}'.format(e))
