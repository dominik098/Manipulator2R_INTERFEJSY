import math
from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

G = 9.8 # przyspieszenie ziemskie
L1 = 1.0 # długość członu pierwszego
L2 = 1.0 # długość członu drugiego
M1 = 3.0 # masa członu pierwszego
M2 = 2.0 # masa członu drugiego
Kp = 1.0 # wzmocnienie członu Kp 

th1_ = 92 # wartość zadana kąta przegubu pierwszego (stopnie)
th2_ = 92 # wartość zadana kąta przegubu drugiego (stopnie)

e1 = 0.0 # wartosc uchybu regulacji przegubu pierwszego
e2 = 0.0 # wartosc uchybu regulacji przegubu drugiego

def derivs(state, t):

    dydx = np. zeros_like(state)

    e1 = th1_ - state[0]
    e2 = th2_ - state[2]
    
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

    dydx[1] = ((+ M22*(Kp*e1 - C12*state[3] - G11) - M12*(Kp*e2 - C21*state[1] - G21))/den)
    
    dydx[2] = state[3]

    dydx[3] = ((- M21*(Kp*e1 - C12*state[3] - G11) + M11*(Kp*e2 - C21*state[1] - G21))/den)
    
    return dydx


# tworzenie tablicy czasu od 0..100 próbkowanej co 0.05 sekundy
dt = 0.05
t = np.arange(0, 20, dt)

# th1 i th2 są początkowymi kątami (stopnie)
# w1 i w2 są początkowymi prędkościami kątowymi (stopnie na sekundę)
th1 = 90
w1 = 0.0
th2 = 90
w2 = 0.0

# stan początkowy
state = np.radians([th1, w1, th2, w2])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, state, t)

x1 = L1*sin(y[:, 0])
y1 = -L1*cos(y[:, 0])

x2 = L2*sin(y[:, 2]) + x1
y2 = -L2*cos(y[:, 2]) + y1

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text


ani = animation.FuncAnimation(fig, animate, range(1, len(y)),
                              interval=dt*1000, blit=True, init_func=init)
plt.show()
