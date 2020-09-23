"""
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 7 - #701 

Newman, Exercise 5.20: A more advanced trapezoidal rule
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

# Part a

def f(x):
    if x == 0: # this would limit the value of the integrand
        return 1
    else:
        return (np.sin(x)/x)**2

def integral_of_f(x):
    return (2 * x * sp.sici(2 * x)[0] + np.cos(2 * x) - 1) / (2 * x)

def g(x):
    return np.log(10 * x + 1)

integration_slices = [] # list of integration slices as an empty array

def integral(f, a, b, error):
    '''
    Adaptive Trapezoidal Method
    Calculate an initial guess (I1), which is one slice.Then, we will double 
    the number of steps, calculating another guess (I2) which is two slices by 
    using the previous guess. Lastly, we will calculate the error between 
    the two guesses. If the error meets the target value, the function returns I2.
    
    error = 1/3 * (I2 - I1)
    target value = (x2 - x1) * delta
    '''
    
    delta = error / (b - a)  # target accuracy per unit interval
    
    def step(x1, x2, f_1, f_2):
        '''
        Parameters:
            
        x1: lower limit of integral
        x2: Upper limit of integral
        f_1 = f(x1): Function evaluated at lower limit
        f_2 = f(x2): Function evaluated at upper limit
        '''
        
        h = np.abs(x2 - x1) # width
        m = np.abs((x1 + x2)/2) # midpoint
        fm = f(m) # integrand at midpoint
        
        I1 = h/2 * (f_1 + f_2) # Trapezoidal Rule
        I2 = h/4 * (f_1 + 2*fm + f_2) # Trapezoidal Rule Estimate
        # I2 = I1/2 + h/2 * fm
        error = np.abs((I2 - I1)/3) 
        S = h/6 * (f_1 + 4*fm + f_2) # Improved Simpson's Rule Estimate
        
        if error < h * delta:
            integration_slices.append(x1)  #adds onto the list
            integration_slices.append(x2)
            return S # returns a value that is more accurate by two orders in h
        else:
            return step(x1, m, f_1, fm) + step(m, x2, fm, f_2) 

    return step(a, b, f(a), f(b))

I = integral(f, 0, 10, 1e-4) 
I_exact = integral_of_f(10)
# I_exact = sp.sici(20)[0] - 1/10 * np.sin(10)**2
error = np.abs(I_exact - I)
print("Part a:")
print('I_exact = {:.16f}'.format(I_exact))
print('I = {:.16f}'.format(I))
print('error = {:.16e}'.format(error))

vals = list(map(f, integration_slices))

##############################################################################

# Part b
'''
Why does the function step(x1 ,x2,f1 ,f2) take not only the positions x1 and x2
as arguments, but also the values f(x1) and f(x2)? Since we know the function
f(x), we could just calculate these values from x1 and x2. Nonetheless, it is a
smart move to include the values of f(x1) and f(x2) as arguments to the function.
Why?
'''
partb = """
Part b: 
If the error does not meet the target, then the function calls itself twice by summing up individual integrals on the 1st and 2nd halves of the interval. 
Doing this will recursively subdivide the integral until it's met by an accurate result. 
Therefore, it is a smart move to include the values f(x1) and f(x2) in the function step(x1, x2, f_1, f_2).
"""

print(partb)

##############################################################################

# Part c
'''
Generate a plot of the integrand with dots added showing where the ends of each 
integration slice lie. You should see larger slices in portions of the integrand
that follow reasonably straight lines (because the trapezoidal rule gives an
accurate value for straight-line integrands) and smaller slices in portions
with more curvature.
'''

textsize = 12
marksize = 3
plt.rc('xtick', labelsize = textsize) 
plt.rc('ytick', labelsize = textsize)

plt.rc('mathtext', fontset = 'stix')
plt.plot(integration_slices, vals, 'o', markersize = marksize, color = "mediumspringgreen")
plt.title("Part C: Modified Plot")
plt.xlabel('x')
plt.ylabel('f(x)')
plt.show()