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
sun = sphere( radius = 2, color = vec(1,1,0) )

programSpeed = 0.00001 
mercury = sphere( radius = mercuryRadius/69173, color = vec(0.9,0.5,1) )
venus = sphere( radius = venusRadius/69173, color = vec(0.8,1,0.5) )
earth = sphere( radius = earthRadius/69173, color = vec(0.5,0.7,1) )
mars = sphere( radius = marsRadius/69173, color = vec(1,0.1,0.8) )
jupiter = sphere( radius = jupiterRadius/69173, color = vec(1,0.8,0.3) )
saturn = sphere( radius = saturnRadius/69173, color = vec(0.9,0.9,1) )

mercuryAngle -= mercuryOrbitRate * programSpeed
venusAngle -= venusOrbitRate * programSpeed
earthAngle -= earthOrbitRate * programSpeed
marsAngle -= marsOrbitRate * programSpeed
jupiterAngle -= jupiterOrbitRate * programSpeed
saturnAngle -= saturnOrbitRate * programSpeed
    
while (True):
    rate(1000)
    
    # change position of each planet
    mercury.pos = vec( mercuryOrbitRadius/155 * cos(radians(mercuryAngle)), 0, mercuryOrbitRadius/155 * sin(radians(mercuryAngle)) ) # using circle trig
    venus.pos = vec( venusOrbitRadius/155 * cos(radians(venusAngle)), 0, venusOrbitRadius/155 * sin(radians(venusAngle)) ) # using circle trig
    earth.pos = vec( earthOrbitRadius/155 * cos(radians(earthAngle)), 0, earthOrbitRadius/155 * sin(radians(earthAngle)) ) # using circle trig
    mars.pos = vec( marsOrbitRadius/155 * cos(radians(marsAngle)), 0, marsOrbitRadius/155 * sin(radians(marsAngle)) ) # using circle trig
    jupiter.pos = vec( jupiterOrbitRadius/155 * cos(radians(jupiterAngle)), 0, jupiterOrbitRadius/155 * sin(radians(jupiterAngle)) ) # using circle trig
    saturn.pos = vec( saturnOrbitRadius/155 * cos(radians(saturnAngle)), 0, saturnOrbitRadius/155 * sin(radians(saturnAngle)) ) # using circle trig
    
    earth.rotate (angle = radians(programSpeed * 360), axis = vec(0, 1, 0)) # rotate the earth 360 times per year
    sun.rotate (angle = radians(programSpeed * 16), axis = vec(0, 1, 0)) # rotate the Sun with a period of about 22 days