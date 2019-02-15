import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

b = 0.01
mu = 0.002
h = 1
tst = 0.01

S0 = 0.999
I0 = 0.001
R0 = 0
b2= 0.01

C = I0 + S0 - 1
l = b - mu + b*C


def xctS(t):
    return 1 + C*(1-mu*t) - l/(b + l*((l-I0*b)/(l*I0*np.exp(b*C/mu)))*np.exp(-l*t + b*C/mu));

def xctI(t):
    return l/(b + l*((l-I0*b)/(l*I0*np.exp(b*C/mu)))*np.exp(-l*t + b*C/mu));

def fS(I, S):
    return -b2*I*S - mu*S + mu

def fI(I, S):
    return b2*I*S -mu*I


Sx = []
Ix = []
T = []
St = []
It = []

S = 0.999
I = 0.001
t = np.linspace(0, 10000, 10000)

for n in range(0, len(t)):

    Sx.append(xctS(t[n]))
    Ix.append(xctI(t[n]))

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


print(len(t))
print(len(T))
print(len(Sx))
print(len(Ix))
print(len(St))
print(len(It))


"""plt.figure(1)
plt.plot(t, S, t, I)
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Exact SIR Solution')
plt.show()"""

plt.plot(t, Sx, 'salmon', t, Ix, 'lawngreen', t, St, 'darkgrey', It, 'darkslategray')
plt.xlabel('Time Step')
plt.ylabel('Population')
plt.title('SIR: Analytic Solution vs Runge-Kutta Method')
plt.legend(labels=['S - Analytic', 'I - Analytic', 'S - Runge-Kutta', 'I - Runge-Kutta'], loc = 'right')
plt.grid()
plt.xlim([0, 2500])
plt.show()







