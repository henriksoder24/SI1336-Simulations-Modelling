import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 12})

#CONSTANTS
N = 1000    #TOTAL POPULATION
nsteps = 10000
b = 0.01
g = 0.001

def fS(I, S):
    return -b*(I*S/N)

def fI(I, S):
    return b*(I*S/N) - g*I

def fR(I):
    return g*I



#INITIAL CONDITIONS
S = [ 900, 800, 700, 600, 500, 400, 300, 200, 100]
I = [100, 200, 300, 400, 500, 600, 700, 800, 900]
R = [ 0, 0, 0, 0, 0, 0, 0, 0, 0]


h = 0.01  #STEP SIZE
t = 0
tst = 0.01

###AMOUNT OF LISTS per ST, IT and RT needed
St = []
It = []
Rt = []
T = []

NlineX = [1000, 0]
NlineY = [0, 1000]

for p in range(0, len(S)):

    St.append([])
    It.append([])
    Rt.append([])
    T.append([])

for n in range(0, len(S)-1):

    A = S[n]
    B = I[n]
    C = R[n]

    for step in range(0, nsteps):

        k1S = fS(B, A)
        k2S = fS(B, A + k1S/2)
        k3S = fS(B, A + k2S/2)
        k4S = fS(B, A + k3S)

        k1I = fI(B, A)
        k2I = fI(B + k1I/2, A)
        k3I = fI(B + k2I/2, A)
        k4I = fI(B + k3I, A)

        k1R = fR(B)
        k2R = fR(B+k1R)
        k3R = fR(B+k2R)
        k4R = fR(B+k3R)

        A += (1/6)*(k1S + 2*k2S + 2*k3S + k4S)
        B += (1/6)*(k1I + 2*k2I + 2*k3I + k4I)
        C += (1/6)*(k1R + 2*k2R + 2*k3R + k4R)

        t += tst
        St[n].append(A)
        It[n].append(B)
        Rt[n].append(C)
        T.append(t)


plt.plot(St[0], It[0], '#1f77b4', St[1], It[1],'#1f77b4', St[2], It[2], '#1f77b4', St[3], It[3], '#1f77b4', St[4], It[4], '#1f77b4', St[5], It[5], '#1f77b4', St[6], It[6], '#1f77b4', St[7], It[7], '#1f77b4', St[8], It[8],'#1f77b4', NlineX, NlineY, '-.')
plt.xlabel('Susceptible')
plt.ylabel('Infected')
plt.title('Epidemic S/I Group Phase Space')
plt.ylim([600, 1000])
plt.xlim([0, 400])
plt.grid()
plt.show()

"""plt.figure(1)
plt.plot(T, St, T, It, T, Rt)
plt.xlabel('Time (Days)')
plt.ylabel('Group Population')
plt.xlim([0, 100])
plt.title('Epidemic SIR Model')
plt.legend(labels=['Susceptible', 'Infected', 'Recovered'], loc = 'right')
plt.grid()
plt.show()"""


