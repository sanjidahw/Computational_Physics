'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 2 - #201

Newman, Exercise 3.5: Visualization of the Solar System
'''
import vpython as vp
import numpy as np

# lists of radius of object(km)
# radius of orbit (km)
# period of orbit (days)
# color
mercury = [2440,57.9e6,88.0, vp.color.blue]
venus = [6052,108.2e6,224.7, vp.color.magenta]
earth = [6371,149.6e6,365.3, vp.color.green]
mars = [3386,227.9e6,687.0, vp.color.red]
jupiter = [69173,778.5e6,4331.6, vp.color.orange]
saturn = [57316,1433.4e6,10759.2, vp.color.cyan]
sun = [695500," "," ", vp.color.yellow]
planets = [mercury, venus, earth, mars, jupiter, saturn]

framerate = 30
timescale = 250
nplanets = 6
scale = 1000
scale_sun = 25

vp.canvas (width=800, height=800)

# start off with the sun
vp.sphere (pos=vp.vector(0.0,0.0,0.0), radius = scale_sun*sun[0], color = sun[3])

# now create the planets
planet = np.empty(nplanets, vp.sphere)
for i in range(nplanets):
    x = planets[i][1]
    y = 0.0
        
    planet[i] = vp.sphere(pos=vp.vector(x,y,0), radius = scale * planets[i][0], color = planets[i][3])

#start time
t = 0.0
while True:
    vp.rate(framerate)   
    t += timescale/framerate
    
# now let's get the planets to move
    for i in range(nplanets):
        x = planets[i][1]*np.cos(2*np.pi*t/planets[i][2])
        y = planets[i][1]*np.sin(2*np.pi*t/planets[i][2])
        planet[i].pos = vp.vector(x, y, 0)
        
earth.rotate (angle = np.radians(0.00001  * 360), axis = vp.vector(0, 1, 0)) # rotate the earth 360 times per year
sun.rotate (angle = np.radians(0.00001  * 16), axis = vp.vector(0, 1, 0)) # rotate the Sun with a period of about 22 days
