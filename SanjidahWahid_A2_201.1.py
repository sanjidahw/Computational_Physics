'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 2 - #201

Newman, Exercise 3.5: Visualization of the Solar System
'''

from vpython import *

# set radius
sunRadius = 695500
mercuryRadius = 2440
venusRadius = 6052
earthRadius = 6371
marsRadius = 3386
jupiterRadius = 69173
saturnRadius = 57316

# set Orbit radius
mercuryOrbitRadius = 57.9
venusOrbitRadius = 108.2
earthOrbitRadius = 149.6
marsOrbitRadius = 227.9
jupiterOrbitRadius = 778.5
saturnOrbitRadius = 1433.4

# set Orbit rate
mercuryOrbitRate = 88.0
venusOrbitRate = 224.7
earthOrbitRate = 365.3
marsOrbitRate = 687.0
jupiterOrbitRate = 4331.6
saturnOrbitRate = 10759.2

# set angle
mercuryAngle = 0
venusAngle = 0
earthAngle = 0
marsAngle = 0
jupiterAngle = 0
saturnAngle = 0
astronomicalUnit = 4.0

# create the Sun object
sun = sphere( radius = sunRadius, color = vec(1,1,0) )

programSpeed = 0.00001 
mercury = sphere( radius = mercuryRadius, color = vec(0.9,0.5,1) )
venus = sphere( radius = venusRadius, color = vec(0.8,1,0.5) )
earth = sphere( radius = earthRadius, color = vec(0.5,0.7,1) )
mars = sphere( radius = marsRadius, color = vec(1,0.1,0.8) )
jupiter = sphere( radius = jupiterRadius, color = vec(1,0.8,0.3) )
saturn = sphere( radius = saturnRadius, color = vec(0.9,0.9,1) )

mercuryAngle -= mercuryOrbitRate * programSpeed
venusAngle -= venusOrbitRate * programSpeed
earthAngle -= earthOrbitRate * programSpeed
marsAngle -= marsOrbitRate * programSpeed
jupiterAngle -= jupiterOrbitRate * programSpeed
saturnAngle -= saturnOrbitRate * programSpeed
    
while (True):
    rate(1000)
    
    # change position of each planet
    mercury.pos = vec( mercuryOrbitRadius * cos(radians(mercuryAngle)), 0, mercuryOrbitRadius * sin(radians(mercuryAngle)) ) # using circle trig
    venus.pos = vec( venusOrbitRadius * cos(radians(venusAngle)), 0, venusOrbitRadius * sin(radians(venusAngle)) ) # using circle trig
    earth.pos = vec( earthOrbitRadius * cos(radians(earthAngle)), 0, earthOrbitRadius * sin(radians(earthAngle)) ) # using circle trig
    mars.pos = vec( marsOrbitRadius * cos(radians(marsAngle)), 0, marsOrbitRadius * sin(radians(marsAngle)) ) # using circle trig
    jupiter.pos = vec( jupiterOrbitRadius * cos(radians(jupiterAngle)), 0, jupiterOrbitRadius * sin(radians(jupiterAngle)) ) # using circle trig
    saturn.pos = vec( saturnOrbitRadius * cos(radians(saturnAngle)), 0, saturnOrbitRadius * sin(radians(saturnAngle)) ) # using circle trig
    
    earth.rotate (angle = radians(programSpeed * 360), axis = vec(0, 1, 0)) # rotate the earth 360 times per year
    sun.rotate (angle = radians(programSpeed * 16), axis = vec(0, 1, 0)) # rotate the Sun with a period of about 22 days