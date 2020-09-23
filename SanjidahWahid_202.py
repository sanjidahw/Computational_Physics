'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 2 - #202

Madelung Constant: Convergence

'''
############################# Part a ##########################################
import numpy as np

n = int(input("Enter a value for L: "))

L_n = list(np.linspace(-n,n,((2*n)+1),True,False,int))    

# we can make use of a function to calculate M(n) for some given value n
def madelung_approx(L_i,L_j,L_k):
    M = 0
    for i in L_i:
        for j in L_j:
            for k in L_k:
                if i == j == k == 0:
                    M += 0
                elif (i+j+k) % 2 == 0:
                    M += (1/np.sqrt(i**2+j**2+k**2))
                else:
                    M += ((-1)*(1/np.sqrt(i**2+j**2+k**2)))
    return(M)
    
madelung_approx(L_n,L_n,L_n)
print("The Madelung constant for the inputted value of L = " + str(n) + " is equal to " + str(madelung_approx(L_n,L_n,L_n)))


# Now we can create a list to calculate the n-1 variant of Madelung constant approximate
L_n_min1 = list(np.linspace(-(n-1),(n-1),((2*(n-1))+1),True,False,int))

# we can make use of a function to calculate delta(n)
def min_diff(L_i,L_j,L_k, n):
    diff = 0
    for i in L_i:
        for j in L_j:
            for k in L_k:
                if (np.abs(i) == n or np.abs(j) == n or np.abs(k) == n) and (i + j + k) % 2 == 0:
                    diff += (1/np.sqrt(i**2+j**2+k**2))
                elif (np.abs(i) == n or np.abs(j) == n or np.abs(k) == n) and (i + j + k) % 2 == 1:
                    diff += ((-1)*(1/np.sqrt(i**2 + j**2 + k**2)))
                else:
                    diff += 0
    return(diff)
    
min_nvalue = madelung_approx(L_n_min1,L_n_min1,L_n_min1)
print("The Madelung constant for the inputted value of L = " + str(n-1) + " is equal to " + str(madelung_approx(L_n_min1,L_n_min1,L_n_min1)))
diff = min_diff(L_n,L_n,L_n,n)
print("M(n) - M(n-1) = " + str(min_diff(L_n,L_n,L_n,n)))
print("M(n-1) + delta(n) = " + str(min_nvalue + diff))

equivalence_boolean = round(madelung_approx(L_n,L_n,L_n), 8) == round(min_nvalue + diff, 8)
print("The Boolean value of M(n) = M(n-1) + delta(n) is about " + str(equivalence_boolean))

# finding Madelung constant with partial sums
def partial_sum(n):
    M = 0
    a = range(1,n + 1)
    for i in a:
        L_i = list(np.linspace(-i,i,((2*i) + 1),True,False,int))
        M += min_diff(L_i,L_i,L_i,i)
    return M

print("The Madelung constant with the use of partial sums and the value of L = " + str(n) + " is " + str(partial_sum(n)))

############################# Part b ##########################################

import matplotlib.pyplot as plt

# let's create a list of values of constant approximation from 0 up to m
def a_values_list(m):
    a = range(m + 1)
    M_list = []
    for i in a:
        L_i = list(np.linspace(-i,i,((2*i) + 1),True,False,int))
        M_list.append(madelung_approx(L_i,L_i,L_i))
    return M_list

def b_values_list(m):
    a = range(m + 1)
    M_list = []
    for i in a:
        L_i = list(np.linspace(-i,i,((2*i) + 1),True,False,int))
        M_list.append((madelung_approx(L_i,L_i,L_i)+(madelung_approx(L_i,L_i,L_i)-min_diff(L_i,L_i,L_i,n)))/2)
    return M_list

# plotting M_n for 10 =< n =< 30
x = np.arange(0, n + 1, 1)

plt.plot(x,a_values_list(n), "-g")
plt.plot(x,b_values_list(n),"-r")

plt.title('M(n) vs. N vs. M(n)=(1/2)(M(n-1) + M(n))')
plt.xlim(10,30) # N values 10....30
