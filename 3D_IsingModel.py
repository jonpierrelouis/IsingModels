#
# 2 Dimensional Ising model
#

import random as rnd, numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate 

def update3d(N, spin, kT, E, M):    # 2D Ising model, N x N lattice
    i, j, k, flip = rnd.randint(0, N-1), rnd.randint(0, N-1), rnd.randint(0, N-1), 0
    dE = 2*spin[i][j][k]*( spin[i-1][j][k] + spin[(i+1)%N][j][k]
                      + spin[i][j-1][k] + spin[i][(j+1)%N][k] 
		      +spin[i][j][k-1] + spin[i][j][(k+1)%N] )
    if (dE < 0.0): flip=1           # flip if dE<0, else flip 
    else:                           # according to exp(-dE/kT)
        p = np.exp(-dE/kT)
        if (rnd.random() < p): flip=1
    if (flip == 1):
        E = E + dE
        M = M - 2*spin[i][j][k]
        spin[i][j][k] = -spin[i][j][k]
    return E, M



def initialize(N):  # set initial spins
    p, spin, E, M = .5, [[[1 for i in range(N)] for j in range(N)] for k in range(N)], 0., 0.     # p = order para.
    for i in range(1, N):
        for j in range(1, N):
            for k in range(1,N):
                if (rnd.random() < p): spin[i][j][k]=-1
                E = E - spin[i-1][j-1][k-1]*spin[i][j][k]
                M = M + spin[i][j][k]

    return spin, E - spin[N-1][N-1][N-1]*spin[0][0][0], M+spin[0][0][0]


#
# MAIN
#

N, passes = 10,10
#N, passes = 1000, 100
iter, Nmc = passes*N**3, passes*N**3
T, Eavg, Mavg = [], [], []
enth = []

for i in range(1,71):               # temperature loop                    
   
    kT = 0.1*i                      # kT = reservoir temperature
    
    spin, E, M  = initialize(N)		
    print(E, M)
    for k in range(iter):           # let it equilibrate
        E, M = update3d(N, spin, kT, E, M)
        
    
    E1, M1 = 0., 0.
    for k in range(Nmc):            # take averages
        E, M = update3d(N, spin, kT, E, M)
        E1, M1 = E1 + E, M1 + M
    E1, M1 = E1/Nmc, M1/Nmc
    T.append(kT), Eavg.append(E1/N**3), Mavg.append(M1/N**3)
    sigma = -1*integrate.trapz(Eavg,T)
    enth.append(sigma)
    #print(sigma, T[-1])

#
#Plot Points
#

plt.figure()
plt.title("Energy")
plt.plot(T, Eavg, 'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle E \rangle/N\epsilon$')
plt.figure()
plt.title("Magnetization")
plt.plot(T, Mavg,'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle M \rangle/N$')

plt.figure()
plt.title("Entropy")
plt.plot(T, enth)
plt.xlabel('$kT/\epsilon$'), plt.ylabel('J/K')

plt.show()


