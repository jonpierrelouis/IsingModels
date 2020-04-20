import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate 

data = np.loadtxt("/home/jpl/testdata_2d.dat")
T = data[:,0]
Eavg = data[:,1]
Mavg = data[:,2]
T_c = 4.4
enth, hc = [], []
for i in range(1,71):
    sigma = -1*integrate.trapz(Eavg[:i],T[:i])
    enth.append(sigma)
    t_i = Eavg[i] - Eavg[i-1]
    t_i = t_i/0.1
    hc.append(t_i)


plt.figure()
plt.plot(T, Eavg, 'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle E \rangle/N\epsilon$')

plt.figure()
plt.plot(T, Mavg,'o')
plt.xlabel('$kT/\epsilon$'), plt.ylabel(r'$\langle M \rangle/N$')

plt.figure()
plt.plot(T, enth)
plt.xlabel('$kT/\epsilon$'), plt.ylabel('S/Nk')

plt.figure()
plt.scatter(T,hc)
plt.xlabel('$kT/\epsilon$'), plt.ylabel('C/Nk')

plt.show()

