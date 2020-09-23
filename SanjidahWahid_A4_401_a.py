'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 4 - #401 

Numerical Convergence of an Integral in Our Solar System
'''

#part (a)
import numpy as np
import matplotlib.pyplot as plt

def f(s):
    return ((1-s)/s)**(1/2)

N_list = [100,200,400,800,1600,3200,6400,12800]
N_array = np.array([100,200,400,800,1600,3200,6400,12800],dtype=int)
N_array = 1/N_array

a = 0.5
b = 1.


I_exact = (np.pi - 2)/4

#Trapezoidal rule evalation for values in list
def trap_eval(lst,fun):
    I_trap_end = 1/2*(fun(a)+fun(b))
    I_trap_val = list(np.zeros(len(lst)))
    for i in range(len(lst)):
        j = int(i)
        N = lst[j]
        h = (b-a)/N
        for n in range(1,N):
            I_trap_val[j] += fun(a+n*h)
        I_trap_val[j] += I_trap_end
        I_trap_val[j] *= h
    return I_trap_val

print(trap_eval(N_list,f))
I_trap_diff = np.array(trap_eval(N_list,f))
I_trap_diff -= I_exact
I_trap_diff = np.abs(I_trap_diff)
print(I_trap_diff)

#Simpson's rule evaluation
def simp_eval(lst,fun):
    I_simp_end = 1/3 * (f(a) + f(b))
    I_simp_val = list(np.zeros(len(lst)))
    for i in range(len(lst)):
        j = int(i)
        N = lst[j]
        h = (b-a)/N
        for n in range(1,N):
            if n%2 == 0:
                I_simp_val[j] += 2/3 * fun(a+n*h)
            else:
                I_simp_val[j] += 4/3 * fun(a+n*h)
        I_simp_val[i] += I_simp_end
        I_simp_val[i] *= h
    return I_simp_val

print(simp_eval(N_list,f))
I_simp_diff = np.array(simp_eval(N_list,f))
I_simp_diff -= I_exact
I_simp_diff = np.abs(I_simp_diff)
print(I_simp_diff)


print(str(I_exact) + " this is exact")

#Plotting for part (a)
fig, (plot1, plot2) = plt.subplots(1,2, figsize=(12,6))
plot1.plot(N_array,I_trap_diff, label = "Trapezoidal rule difference")
plot1.plot(N_array,I_simp_diff, label = "Simpson's rule difference")
plot1.set_title(r"Log scale of $\Delta$ versus $\frac{1}{N}$ for Part (a)")
plot1.set_ylabel(r"$\frac{1}{N}$", fontsize = 14)
plot1.set_xlabel(r"$\Delta$", fontsize = 14)
plot1.loglog()
plot1.legend()

##############################################################################
#part (b)

const = (b**(-.5)*(1-b)**(3/2)) - (a**(-.5)*(1-a)**(3/2))

def alt_f(s):
    return (((1-s)/s)**(3/2))

#Trapezoidal rule evalation for values in list, modified for alt_f
def trap_eval_alt(lst):
    I_trap_end = 1/2*(alt_f(a)+alt_f(b))
    I_trap_val = list(np.zeros(len(lst)))
    for i in range(len(lst)):
        j = int(i)
        N = lst[j]
        h = (b-a)/N
        for n in range(1,N):
            I_trap_val[j] += alt_f(a+n*h)
        I_trap_val[j] += I_trap_end
        I_trap_val[j] *= h/2
        I_trap_val[j] += const
        I_trap_val[j] *= -2/3
    return I_trap_val

I_trap_diff_alt = np.array(trap_eval_alt(N_list))
I_trap_diff_alt -= I_exact
I_trap_diff_alt = np.abs(I_trap_diff_alt)
print(I_trap_diff_alt)

#Simpson's rule evaluation, modified for alt_f
def simp_eval_alt(lst):
    I_simp_end = 1/3 * (alt_f(a) + alt_f(b))
    I_simp_val = list(np.zeros(len(lst)))
    for i in range(len(lst)):
        j = int(i)
        N = lst[j]
        h = (b-a)/N
        for n in range(1,N):
            if n%2 == 0:
                I_simp_val[j] += 2/3 * alt_f(a+n*h)
            else:
                I_simp_val[j] += 4/3 * alt_f(a+n*h)
        I_simp_val[i] += I_simp_end
        I_simp_val[i] *= h/2
        I_simp_val[i] += const
        I_simp_val[i] *= -2/3
    return I_simp_val

I_simp_diff_alt = np.array(simp_eval_alt(N_list))
I_simp_diff_alt -= I_exact
I_simp_diff_alt = np.abs(I_simp_diff_alt)
print(I_simp_diff_alt)

#Plotting for part (b)
plot2.plot(N_array,I_trap_diff_alt, label = "Trapezoidal rule difference")
plot2.plot(N_array,I_simp_diff_alt, label = "Simpson's rule difference")
plot2.set_title(r"Log scale of $\Delta$ versus $\frac{1}{N}$ for Part (b)")
plot2.set_ylabel(r"$\frac{1}{N}$", fontsize = 14)
plot2.set_xlabel(r"$\Delta$", fontsize = 14)
plot2.loglog()
plot2.legend()
fig.show()

##############################################################################
#part (c)




