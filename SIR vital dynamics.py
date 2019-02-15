import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

#CONSTANTS
N = 1000    #TOTAL POPULATION
nsteps = 20000
b = 0.01
g = 0.001
mu = 0.00005

def fS(I, S):
    return mu*N - mu*S -b*(I*S/N)

def fI(I, S):
    return b*(I*S/N) - g*I - mu*I

def fR(I):
    return g*I - mu*R

St = [ ]
It = [ ]
Rt = [ ]
T = [ ]

#INITIAL CONDITIONS
S = 999
I = 1
R = 0

h = 0.01  #STEP SIZE
t = 0
tst = 0.01

for step in range(0, nsteps):

    k1S = fS(I, S)
    k2S = fS(I, S + k1S/2)
    k3S = fS(I, S + k2S/2)
    k4S = fS(I, S + k3S)

    k1I = fI(I, S)
    k2I = fI(I + k1I/2, S)
    k3I = fI(I + k2I/2, S)
    k4I = fI(I + k3I, S)

    k1R = fR(I)
    k2R = fR(I+k1R)
    k3R = fR(I+k2R)
    k4R = fR(I+k3R)

    S += (1/6)*(k1S + 2*k2S + 2*k3S + k4S)
    I += (1/6)*(k1I + 2*k2I + 2*k3I + k4I)
    R += (1/6)*(k1R + 2*k2R + 2*k3R + k4R)

    t += tst
    St.append(S)
    It.append(I)
    Rt.append(R)
    T.append(t)

plt.figure(1)
plt.plot(T, St, T, It, T, Rt)
plt.xlabel('Time (Days)')
plt.ylabel('Group Population')
#plt.xlim([0, 100])
plt.title('Endemic SIR Model')
plt.legend(labels=['Susceptible', 'Infected', 'Recovered'], loc = 'right')
plt.grid()
plt.show()

plt.figure(2)
plt.plot(St, It)
plt.xlabel('Susceptible')
plt.ylabel('Infected')
plt.title('Susceptible/Infected Group Phase Space')
plt.grid()
plt.show()
