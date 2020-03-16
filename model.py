import math
from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

G = 9.8 # przyspieszenie ziemskie
L1 = 1.0 # długość członu pierwszego
L2 = 1.0 # długość członu drugiego
M1 = 1.0 # masa członu pierwszego
M2 = 1.0 # masa członu drugiego
Kp = 1.0 # wzmocnienie członu Kp 

def derives(state, t):

    dydx = np. zeros_like(state)

    suma = state[0] + state[2]
    
    M11 = (M1 + M2)*L1**2 + M2*L2**2 + 2*M2*L1*L2*cos(state[2])
    M12 = M2*L2**2 + M2*L1*L2*cos(state[2])
    M21 = M2*L2**2 + M2*L1*L2*cos(state[2])
    M22 = M2*L2**2

    C12 = -2*M2*L1*L2*sin(state[2])*state[1]*state[3] - M2*L1*L2*sin(state[2])*state[3]*state[3]
    C21 = M2*L1*L2*sin(state[2])*state[1]*state[1]

    G11 = (M1+M2)*G*L1*sin(state[0]) + M2*G*L2*sin(suma)
    G21 = M2*G*L2*sin(suma)

    den = M11*M22 - M12*M21
    
    dydx[0] = state[1]

    dydx[1] = ( + M22*(Kp*state[0] - C12*state[3]
                - G11) - M12*(Kp*state[2]
                - C21*state[1] - G21))/den)
    
    dydx[2] = state[3]

    dydx[3] = ( - M21*(Kp*state[0] - C12*state[3]
                - G11) + M11(Kp*state[2]
                - C21*state[1] - G21)/den)
    
    return dydx


# tworzenie tablicy czasu od 0..100 próbkowanej co 0.05 sekundy
dt = 0.05
t = np.arrange(0, 20, dt)

# th1 i th2 są początkowymi kątami (stopnie)
# w1 i w2 są początkowymi prędkościami kątowymi (stopnie na sekundę)
th1 = 0.0
w1 = 0.0
th2 = 0.0
w2 = 0.0

# stan początkowy
state = np.radians([th1, w1, th2, w2])
