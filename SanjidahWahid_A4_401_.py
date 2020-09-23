'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 4 - #401 

Numerical Convergence of an Integral in Our Solar System
'''

# Part a

import numpy as np
import matplotlib.pyplot as plt


# Trapezoidal Rule

a = 1/2
b = 1
M = 100
N_list = [100,200,400,800,1600,3200,6400,12800]

def g(s):
	return np.sqrt((1 - s)/s)

def I_T(M):
    I_T_end = 1/2 * (g(a) + g(b))
    h = (b-a)/M
    t_odd = 0
    t_even = 0
    for n in range (1, M//2):
        t_odd += g(a + (2*n - 1)*h)
        t_even += g(a + 2*n*h)
    t_odd += g(b - h)
    t_int = h*(I_T_end + t_odd + t_even)
    return t_int

print("Part a: With the use of Trapezoidal Rule and with N =  " + str(M) + ", the numerical evaluation is " + str(I_T(M)))

#Simpson's Rule

def I_S(M):
    I_S_end = 1/3 *(g(a) + g(b))
    h = (b - a)/M
    s_odd = 0
    s_even = 0
    for n in range (1, M//2):
        s_odd += g(a + (2*n - 1)*h)
        s_even += g(a + 2*n*h)
    s_odd += g(b - h)
    s_int = h*(I_S_end + (4/3) * s_odd + (2/3) * s_even)
    return s_int

print("Part a: With the use of Simpson's Rule and with N =  " + str(M) + ", the numerical evaluation is " + str(I_S(M)))


Trap_list = []
Simp_list = []
N_array = np.array(N_list)
N_arraya1 = 1/N_array
I_exact = (np.pi - 2)/4


for i in range(len(N_list)):
    Simp_list.append(np.abs(I_S(N_list[i])-I_exact))
    print(I_S(N_list[i]))
    
for i in range(len(N_list)):
    Trap_list.append(np.abs(I_T(N_list[i])-I_exact))
    print(I_T(N_list[i]))
	
fig_a = plt.figure(figsize = (12,6))
plt.plot(N_arraya1,Trap_list, label = "Part a: Trapezoidal rule absolute difference", color = "darkgreen")
plt.plot(N_arraya1,Simp_list, label = "Part a: Simpson's rule absolute difference", color = "firebrick")
plt.title(r"Part a: Log scale of $\Delta$ versus $\frac{1}{N}$")
plt.ylabel(r"$\frac{1}{N}$", fontsize = 12)
plt.xlabel(r"$\Delta$", fontsize = 12)
plt.loglog()
plt.legend()

###############################################################################

# Part b

# Trapezoidal Rule

a = 1/2
b = 1

def gb(s):
	return ((1 - s)/s)**(3/2)

def fb(s):
	return (1 - s)**(3/2) /np.sqrt(s)

def I_T_b(Mb):
    I_Tb_end = 1/2 *(gb(a) + gb(b))
    h = (b - a)/Mb
    tb_odd = 0
    tb_even = 0
    for n in range (1, Mb//2):
        tb_odd += gb(a + (2*n - 1)*h)
        tb_even += gb(a + 2*n*h)
    tb_odd += gb(b - h)
    t_b_int = -2/3 *((1/2)*h*(I_Tb_end + tb_odd + tb_even) + (fb(b) - fb(a)))
    return t_b_int

print("Part b: With the use of Trapezoidal Rule and with N =  " + str(M) + ", the numerical evaluation is " + str(I_T_b(M)))

#Simpson's Rule

def I_S_b(Mb):
    I_Sb_end = 1/3 *(gb(a) + gb(b))
    h = (b - a)/Mb
    sb_odd = 0
    sb_even = 0
    for n in range (1, Mb//2):
        sb_odd += gb(a + (2*n - 1)*h)
        sb_even += gb(a + 2*n*h)
    sb_odd += gb(b - h)
    s_b_int = -2/3*((1/2)*h*(I_Sb_end + (4/3) * sb_odd + (2/3) * sb_even)+ (fb(b) - fb(a)))
    return s_b_int

print("Part b: With the use of Simpson's Rule and with N =  " + str(M) + ", the numerical evaluation is " + str(I_S_b(M)))

N_list = [100,200,400,800,1600,3200,6400,12800]
Trap_listb = []
Simp_listb = []
N_array = np.array(N_list)
N_arrayb2 = 1/N_array
I_exact = (np.pi - 2)/4

for i in range(len(N_list)):
    Simp_listb.append(np.abs(I_S_b(N_list[i])-I_exact))
    print(I_S_b(N_list[i]))
    
for i in range(len(N_list)):
    Trap_listb.append(np.abs(I_T_b(N_list[i])-I_exact))
    print(I_T_b(N_list[i]))
	
fig_b = plt.figure(figsize = (12,6))
plt.plot(N_arrayb2,Trap_listb, label = "Part b: Trapezoidal rule absolute difference", color = "steelblue")
plt.plot(N_arrayb2,Simp_listb, label = "Part b: Simpson's rule absolute difference", color = "orchid")
plt.title(r"Part b: Log scale of $\Delta$ versus $\frac{1}{N}$")
plt.ylabel(r"$\frac{1}{N}$", fontsize = 12)
plt.xlabel(r"$\Delta$", fontsize = 12)
plt.loglog()
plt.legend()

###############################################################################

# Part c

#defining the y = ln(delta) of data sets

#Trapezoidal
ln_T_list = np.log(Trap_list)
ln_T_listb = np.log(Trap_listb)

#Simpson's
ln_S_list = np.log(Simp_list)
ln_S_listb = np.log(Simp_listb)

ln_N_array = np.log(N_array)

def coeff(x, y):
	n = len(x)
	mean_x, mean_y = np.mean((x)), np.mean((y))
	xy = np.sum((y)*(x)) - n*mean_y*mean_x
    
	xx = np.sum((x)*(x)) - n*mean_x*mean_x
	alpha = xy/xx
	beta = mean_y - alpha*mean_x
	return alpha, beta

print(str(coeff(ln_N_array, ln_T_list)) + " values for alpha and beta for part a Trapezoidal")
print(str(coeff(ln_N_array, ln_T_listb)) + " values for alpha and beta for part b Trapezoidal")
print(str(coeff(ln_N_array, ln_S_list)) + " values for alpha and beta for part a Simpson's")
print(str(coeff(ln_N_array, ln_S_listb)) + " values for alpha and beta for part b Simpson's")

#linear approximation function
def lin_reg(N, a, b):
	return a*N + b

R_2 = []        #R^2
error = [ln_T_list, ln_S_list, ln_T_listb, ln_S_listb]

for i in range (4):
	corr_m = np.corrcoef(ln_N_array, error[i])
	corr_xy = corr_m[0,1]
	R_2.append(corr_xy**2)

print(str(R_2[0]) + " : Trapezoidal part a R^2")
print(str(R_2[1]) + " : Trapezoidal part b R^2")
print(str(R_2[2]) + " : Simpson's part a R^2")
print(str(R_2[3]) + " : Simpson's part b R^2")	

fig_c = plt.figure(figsize = (12,6))
plt.plot(ln_N_array, ln_S_list, ".", color = "red", label = str(R_2[2]) + " : Simpson's part a R^2, Best Fit = " + str(coeff(ln_N_array, ln_S_list)[0]) + r"$\bullet x + $" + str(coeff(ln_N_array, ln_S_list)[1]))
plt.plot(ln_N_array, lin_reg(ln_N_array, coeff(ln_N_array, ln_S_list)[0], coeff(ln_N_array,ln_S_list)[1]), "--", color = "maroon")
plt.plot(ln_N_array, ln_S_listb, ".", color = "cyan", label = str(R_2[3]) +" : Simpson's part b R^2, Best Fit = " + str(coeff(ln_N_array,ln_S_listb)[0]) + r"$\bullet x + $" + str(coeff(ln_N_array,ln_S_listb)[1]))
plt.plot(ln_N_array, lin_reg(ln_N_array, coeff(ln_N_array, ln_S_listb)[0], coeff(ln_N_array,ln_S_listb)[1]), "--", color = "midnightblue")
plt.plot(ln_N_array, ln_T_list, ".", color = "limegreen", label = str(R_2[0]) + " : Trapezoidal part a R^2, Best Fit = " + str(coeff(ln_N_array,ln_T_list)[0]) + r"$\bullet x + $" + str(coeff(ln_N_array,ln_T_list)[1]))
plt.plot(ln_N_array, lin_reg(ln_N_array, coeff(ln_N_array, ln_T_list)[0], coeff(ln_N_array,ln_T_list)[1]), "--", color = "palevioletred")
plt.plot(ln_N_array, ln_T_listb, ".", color = "gold", label = str(R_2[1]) + " : Trapezoidal part b R^2, Best Fit = " + str(coeff(ln_N_array,ln_T_listb)[0]) + r"$\bullet x + $" + str(coeff(ln_N_array,ln_T_listb)[1]))
plt.plot(ln_N_array, lin_reg(ln_N_array, coeff(ln_N_array, ln_T_listb)[0], coeff(ln_N_array,ln_T_listb)[1]), "--", color = "lightseagreen")
plt.title(r"Part c: Linear Regression for ln$\Delta$ Datapoints")
plt.legend(fontsize = 8)

###############################################################################

# Part d

def gaussxw(N): 

	a = np.linspace(3, 4*N - 1, N)/(4*N + 2) 
	x = np.cos(np.pi*a + 1/(8*N*N*np.tan(a)))
	epsilon = 1e-15 
	delta = 1.0
	
	while delta>epsilon:
		p0 = np.ones(N,float)
		p1 = np.copy(x)
		
		for k in range(1,N):
			p0,p1 = p1,((2*k + 1)*x*p1 - k*p0)/(k + 1) 
		dp = (N + 1)*(p0 - x*p1)/(1 - x*x) 
		dx = p1/dp 
		x -= dx 
		delta = max(abs(dx), default = 0) 
	
	# weights 
	w = 2*(N + 1)*(N + 1)/(N*N*(1 - x*x)*dp*dp) 
	return x,w
 
def gaussxwab(N,a,b): 
	x,w = gaussxw(N) 
	return 0.5*(b - a)*x + 0.5*(b + a),0.5*(b - a)*w

N = np.linspace(25,150,26)

gaussq_list = []
for i in N:
	x,w = gaussxwab(int(i),a,b)
	I = np.sum(w*g(x))
	gaussq_list.append(I)
 
gaussq_array = np.array(gaussq_list)
gaussq_array -= I_exact
gaussq_array = np.abs(gaussq_array)
print(gaussq_array)

N = 1/N

# find coeffs alpha, beta, and R^2 for ln Delta and ln N
ln_gaussq = np.log(gaussq_array) # y values
ln_N = np.log(np.linspace(25,150,26)) # x values


R_2_gaussq = (np.corrcoef(ln_N, ln_gaussq)[0][1])**2

print(str(coeff(ln_N, ln_gaussq)) + " values for alpha and beta for Gaussian Quadrature")

fig_d, (plot3, plot4) = plt.subplots(1,2, figsize=(12,6))
plot3.plot(N, gaussq_array, "x", color = "orange", label = "Gaussian Quadrature from Eqn. 1")
plot3.loglog()
plot3.set_title(r"$\Delta$ versus $\frac{1}{N}$: Part d")
plot3.set_xlabel(r"$\frac{1}{N}$", fontsize = 14)
plot3.set_ylabel(r"$\Delta$", fontsize = 14)
plot3.legend()

fig_lin_reg_gauss = plt.figure(figsize=(12,6))
plt.plot(ln_N, ln_gaussq, "c.", label = str(R_2_gaussq) + " : Gaussian Quadrature R^2, Best Fit = " + str(coeff(ln_N, ln_gaussq)[0]) + r"$\bullet x + $" + str(coeff(ln_N, ln_gaussq)[1]))
plt.plot(ln_N, lin_reg(ln_N, coeff(ln_N, ln_gaussq)[0],coeff(ln_N, ln_gaussq)[1]), "c--")
plt.title(r"Part d: Linear Regression for ln$\Delta$ Datapoints")
plt.legend()


###############################################################################

# Part e

gaussq_e = []
for i in N:
    x,w = gaussxwab(int(i),a,b)
    I = (np.sum(w*gb(x)) * -1/3) + 1/3
    gaussq_e.append(I)
	
gaussq_e = np.array(gaussq_e)
gaussq_e -= I_exact
gaussq_e = np.abs(gaussq_e)
print(gaussq_e)

N = np.linspace(25,150,26)
N = 1/N

plot4.plot(N, gaussq_array,"rx", label = "Gaussian Quadrature from Eqn. 2")
plot4.loglog()
plot4.set_title(r"$\Delta$ versus $\frac{1}{N}$: part e")
plot4.set_xlabel(r"$\frac{1}{N}$", fontsize = 14)
plot4.set_ylabel(r"$\Delta$", fontsize = 14)
plot4.legend()

# find alpha, beta, and R^2 for ln Delta and ln N
ln_gaussq_e = np.log(gaussq_e) # y values
ln_N = ln_N # x values

R_2_gaussq_e = (np.corrcoef(ln_N, ln_gaussq_e)[0][1])**2

fig_gauss_alt = plt.figure(figsize=(12,6))
plt.plot(ln_N, ln_gaussq_e, "m.", label = str(R_2_gaussq_e) + " : Equation 2 Gaussian Quadrature R^2, Best Fit = " + str(coeff(ln_N, ln_gaussq_e)[0]) + r"$\bullet x + $" + str(coeff(ln_N, ln_gaussq_e)[1]))
plt.plot(ln_N, lin_reg(ln_N, coeff(ln_N, ln_gaussq_e)[0], coeff(ln_N, ln_gaussq_e)[1]), "m--")
plt.title(r"Part e: Linear Regression for ln$\Delta$ Datapoints")
plt.legend()