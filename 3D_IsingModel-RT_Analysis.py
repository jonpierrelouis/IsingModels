#
# 2 Dimensional Ising model
#

import random as rnd, numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate 
#from mpl_toolkits import mplot3d # noqa: F401 unused import
#import matplotlib.pyplot as plt
#rimport vpython as vp
#from vpython import *

def update3d(N, spin, kT, E, M, arr):    # 2D Ising model, N x N lattice
    i, j, k, flip = rnd.randint(0, N-1), rnd.randint(0, N-1), rnd.randint(0, N-1), 0
    dE = 2*spin[i][j][k]*( spin[i-1][j][k] + spin[(i+1)%N][j][k]
                      + spin[i][j-1][k] + spin[i][(j+1)%N][k] 
		      +spin[i][j][k-1] + spin[i][j][(k+1)%N] )
    if (dE < 0.0): flip=1           # flip if dE<0, else flip 
    else:                           # according to exp(-dE/kT)
        p = np.exp(-dE/kT)
        #print(i, j, k, arr[i,j,k])
        if (rnd.random() < p): flip=1
    if (flip == 1):
        E = E + dE
        M = M - 2*spin[i][j][k]
        spin[i][j][k] = -spin[i][j][k]
        arr[i,j,k] = spin[i][j][k]
        #print(i, j, k, arr[i,j,k])
    return E, M, arr

def initialize(N):  # set initial spins
    p, spin, E, M = 0.5, [[[1 for i in range(N)] for j in range(N)] for k in range(N)], 0., 0.     # p = order para.
    arr = np.zeros((N,N,N))    
    
    for i in range(1, N):
        for j in range(1, N):
            for k in range(1, N):
                if (rnd.random() < p): spin[i][j][k]=-1
                #arr[i][j][k] = spin[i][j][k]
                E = E - spin[i-1][j-1][k-1]*spin[i][j][k]
                M = M + spin[i][j][k]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                arr[i,j,k] = spin[i][j][k]
                #print(i, j, k, arr[i,j,k])

    return spin, E - spin[N-1][N-1][N-1]*spin[0][0][0], M+spin[0][0][0], arr

N, passes = 100, 100
#N, passes = 1000, 100
iter, Nmc = passes*N**2, passes*N**2
T, Eavg, Mavg = [], [], []
enth = []
canvas(title='3D Spin alignment',
     width=800, height=400,
     center=vector(N/2,N/2,N/2), background=color.black)

#for i in range(1,41):               # temperature loop                    
for i in range(1,50):   
    kT = 0.1*i                      # kT = reservoir temperature
    
    spin, E, M, in_arr  = initialize(N)
    #print(in_arr)
    for k in range(iter*10):           # let it equilibrate
        E, M, spin_arr = update3d(N, spin, kT, E, M,in_arr)  
 
    E1, M1 = 0., 0.
    for k in range(Nmc*10):            # take averages
        E, M, spin_arr = update3d(N, spin, kT, E, M, in_arr)
        E1, M1 = E1 + E, M1 + M
    E1, M1 = E1/Nmc, M1/Nmc
    T.append(kT), Eavg.append(E1/N), Mavg.append(M1/N)
    sigma = -1*integrate.trapz(Eavg,T)
    enth.append(sigma)
    #print(sigma, T[-1])


plt.figure()
plt.plot(T, Eavg, 'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle E \rangle/N\epsilon$')
plt.figure()
plt.plot(T, Mavg,'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle M \rangle/N$')

plt.figure()
plt.plot(T, enth)
plt.xlabel('Temperature'), plt.ylabel('Entropy')
plt.show()

plt.figure()
plt.axes(projection='3d')
plt.scatter(spin_arr[0],spin_arr[1], spin_arr[2])
plt.show()



canvas(title='3D Spin alignment',
     width=800, height=400,
     center=vector(N/2,N/2,N/2), background=color.black)
for i in range(N):
    for j in range(N):
        for k in range(N):
            if (spin_arr[i,j,k] == 1):
                sphere(pos=vector(i,j,k), radius=0.25, color=color.yellow)
            else:
                sphere(pos=vector(i,j,k), radius=0.25, color=color.red)


