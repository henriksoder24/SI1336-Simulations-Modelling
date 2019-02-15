import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

#CONSTANTS
N = 1000    #TOTAL POPULATION
nsteps = 10000
b = 0.01
g = 0.0005

def fS(I, S):
    return -b*(I*S/N)

def fI(I, S):
    return b*(I*S/N) - g*I


St = [ ]
It = [ ]

T = [ ]

#INITIAL CONDITIONS
S = 999
I = 1


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


    S += (1/6)*(k1S + 2*k2S + 2*k3S + k4S)
    I += (1/6)*(k1I + 2*k2I + 2*k3I + k4I)

    t += tst
    St.append(S)
    It.append(I)
    T.append(t)

plt.figure(1)
plt.plot(T, St, T, It)
plt.xlabel('Time (Days)')
plt.ylabel('Group Population')
plt.xlim([0, 100])
plt.title('Epidemic SIR Model')
plt.legend(labels=['Susceptible', 'Infected', 'Recovered'], loc = 'right')
plt.grid()
plt.show()


##Statistical Analysis: VARIANCE
print(np.var(St))
print(np.var(It))


"""plt.figure(2)
plt.plot(St, It)
plt.xlabel('Susceptible')
plt.ylabel('Infected')
plt.title('Susceptible/Infected Group Phase Space')
plt.grid()
plt.show()"""
