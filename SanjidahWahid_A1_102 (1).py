'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 1 - #102

Newman, Exercise 2.9: Madelung Constant with the following modification:
    For ease of comparison, forego the issue of convergence and choose L= 25a.  
    Have your program print the answer, which is M=−1.77020578···.
'''

from numpy import sqrt

L = 25
M = 0

## slow

for i in range(-L, L+1):
    for j in range(-L,L+1):
        for k in range(-L, L+1):
            if i != 0 or j != 0 or k != 0:
                M += ((-1)**(i+j+k))/sqrt(i**2+j**2+k**2)
                
print(f'The value of the Madelung constant is {M} with {L} atoms in each direction.')