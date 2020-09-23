"""
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 6 - #601 

Newman, Exercise 8.10: Cometary Orbits
"""

from numpy import sqrt, array, empty
#import matplotlib.pyplot as plt

# Parameters
G = 6.67e-11
M = 1.9891e30 # Mass of the sun (kg)
t0 = 0.0 # intial time (s)
tf = 2e9 # final time (s)
H = tf # Initial step size
n_max = 8 # Maximum allowed # of splits of the current (sub)interval of integration
delta =  1000 / (3600*24*365.25) # accuracy in position per unit of time (m/s)

# Initial conditions
x0 = 4e12 # m
y0 = 0.0
vx0 = 0.0
vy0 = 500.0  # (m/s)
r0 = array([x0, vx0, y0, vy0], float) 

# Returns the components of the current acceleration and velocity vectors
def f(r):
    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]
    dist = sqrt(x ** 2 + y ** 2)
    ax = -G * M * x / dist ** 3
    ay = -G * M * y / dist ** 3
    return array([vx, ax, vy, ay], float)

# Lists for storing positions at different time t
xpoints = [x0]
ypoints = [y0]
tpoints = [t0]

def step(r,t,H):
    n = 1 
    r1 = r + 0.5*H*f(r)
    r2 = r + H*f(r)

    # Richardson extrapolation
    Rnew = empty([1,4],float)
    Rnew[0] = 0.5*(r1 + r2 + 0.5*H*f(r2))

    error = 2*H*delta
    while error > H*delta:
        if n < n_max:
            n += 1
            h = H/n
            r1 = r + 0.5*h*f(r)
            r2 = r + h*f(r)
            for i in range(n - 1):
                r1 += h*f(r2)
                r2 += h*f(r1)

            Rold  = Rnew 
            Rnew  = empty([n,4], float)
            Rnew[0] = 0.5*(r1 + r2 + 0.5 * h * f(r2))
            for m in range(1,n):
                # Error estimate to O(h^(2m))
                epsilon = (Rnew[m-1] - Rold[m-1]) / ((n / (n - 1))**(2*m) - 1)
                Rnew[m] = Rnew[m-1] + epsilon
            # Euclidian error estimate
            error = sqrt(epsilon[0] ** 2 + epsilon[2] ** 2)
       
        else: break

    if error > H * delta and n > n_max:
        R1 = step(r, t, H / 2)
        R2 = step(r, t + H/2, H/2)
    else:
        r += Rnew[n-1]
        xpoints.append(Rnew[n-1, 0]) # Best estimate of x at t+h
        ypoints.append(Rnew[n-1, 2]) # Best estimate of y at t+h
        tpoints.append(t+h)

