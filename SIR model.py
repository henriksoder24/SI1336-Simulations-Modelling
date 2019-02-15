import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

#CONSTANTS
N = 1000    #TOTAL POPULATION
nsteps = 10000
b = 0.010005
g = 0.00050005

def fS(I, S):
    return -b*(I*S/N)

def fI(I, S):
    return b*(I*S/N) - g*I

def fR(I):
    return g*I

St = [ ]
It = [ ]
Rt = [ ]
T = [ ]

#INITIAL CONDITIONS
S = 996
I = 1
R = 3

h = 1  #STEP SIZE
t = 0
tst = 0.01

for step in range(0, nsteps):

    k1S = fS(I, S)*h
    k2S = fS(I, S + k1S/2)*h
    k3S = fS(I, S + k2S/2)*h
    k4S = fS(I, S + k3S)*h

    k1I = fI(I, S)*h
    k2I = fI(I + k1I/2, S)*h
    k3I = fI(I + k2I/2, S)*h
    k4I = fI(I + k3I, S)*h

    k1R = fR(I)*h
    k2R = fR(I+k1R)*h
    k3R = fR(I+k2R)*h
    k4R = fR(I+k3R)*h

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
plt.xlim([0, 100])
plt.title('Epidemic SIR Model')
plt.legend(labels=['Susceptible', 'Infected', 'Recovered'], loc = 'right')
plt.grid()
plt.show()

"""plt.figure(2)
plt.plot(St, It)
plt.xlabel('Susceptible')
plt.ylabel('Infected')
plt.title('Susceptible/Infected Group Phase Space')
plt.grid()
plt.show()"""

print(np.var(St))
print(np.var(It))
print(np.var(Rt))

print(St)
