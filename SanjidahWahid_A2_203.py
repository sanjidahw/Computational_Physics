'''
Sanjidah Wahid
PHYS 39906: Computational Physics Summer 2020
Assignment 2 - #203

Newman, Exercise 3.7: Mandelbrot Set with modifications

'''

import numpy as np
import matplotlib.pyplot as plt

xmax = 2 
N = 100

x=np.linspace(-xmax, xmax, N)
y=x
xv, yv = np.meshgrid(x, y)
mandelbrot=np.zeros([N,N])
niter = 100
c = 1


for i in range(N):
    for j in range(N):
        z = 0
        for k in range(niter):
            z = z**2+xv[i,j]+1j*yv[i,j]
            if np.abs(z)>2:
                mandelbrot[i,j] = 1
                break
            else:
                mandelbrot[i,j] = 0

plt.subplot(121)
plt.title('Mandelbrot set')
im=plt.imshow(mandelbrot, origin='lower', extent=[-xmax, xmax, -xmax, xmax], cmap='gray')
plt.xlabel('Real part of c')
plt.ylabel('Imaginary part of c')
clb=plt.colorbar(im,fraction = 0.046, pad = 0.04)
clb.ax.set_ylabel('# of iterations', rotation = 270)

mandelbrot=np.zeros([N,N])
for i in range(N):
    for j in range(N):
        z = 0
        for k in range(niter):
            z = z**2+xv[i,j]+1j*yv[i,j]
            if np.abs(z)>2:
                mandelbrot[i,j] = k
                break
            else:
                mandelbrot[i,j] = niter

plt.subplot(122)
plt.title('Mandelbrot image of convergence')
im=plt.imshow(mandelbrot, origin='lower', extent=[-xmax, xmax, -xmax, xmax], cmap='jet')
plt.xlabel('Real part of c')
plt.ylabel('Imaginary part of c')
clb=plt.colorbar(im,fraction = 0.046, pad = 0.04)
clb.ax.set_ylabel('# of iterations', rotation = 270)

plt.tight_layout()
