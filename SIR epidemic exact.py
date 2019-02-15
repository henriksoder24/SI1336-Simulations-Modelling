import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

b = 0.01
mu = 0.001

"""S = 0.0
I = 0.0
R = 0.0"""

S0 = 999
I0 = 1
R0 = 0

C = I0 + S0 - 1
l = b - mu + b*C

def xctS(t):
    return 1 + C*(1-mu*t) - l/(b + l*((l-I0*b)/(l*I0*np.exp(b*C/mu)))*np.exp(-l*t + b*C/mu))



def xctI(t):
    return l/(b + l*((l-I0*b)/(l*I0*np.exp(b*C/mu)))*np.exp(-l*t + b*C/mu))

S = []
I = []

t = np.linspace(0, 100, 10000)

for n in range(0, len(t)):

    S.append(xctS(t[n]))
    I.append(xctI(t[n]))

plt.figure(1)
plt.plot(t, S, t, I)
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('Exact SIR Solution')
plt.show()







