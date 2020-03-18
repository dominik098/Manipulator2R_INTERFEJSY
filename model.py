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

Kp1 = 0.35 # wzmocnienie członu proporcjonalnego przegub 1
Ki1 = 0.05 # wzmocnienie członu całkującego przegub 1
Kd1 = 0.008 # wzmocnienie członu różniczukjącego przegub 1

Kp2 = 0.38 # wzmocnienie członu proporcjonalnego przegub 2
Ki2 = 0.2 # wzmocnienie członu całkującego przegub 2
Kd2 = 0.0001 # wzmocnienie członu różniczukjącego przegub 2

th1_ = 30.0 # wartość zadana kąta przegubu pierwszego (stopnie)
th2_ = 10.0 # wartość zadana kąta przegubu drugiego (stopnie)

th1_r = np.radians(th1_)
th2_r = np.radians(th2_)

e1 = 0.0 # wartosc uchybu regulacji przegubu pierwszego
e2 = 0.0 # wartosc uchybu regulacji przegubu drugiego

uchyb_poprzedni1 = 0.0 # zmienna buforu uchybu pierwszego
uchyb_poprzedni2 = 0.0 # zmienna buforu uchyba drugiego

def derivs(state, t):

    dydx = np. zeros_like(state)
    
    # regulator PID przegubu pierwszego
    e1 = th1_ - state[0]
    pochodna1 = (e1 - uchyb_poprzedni1)/dt
    sterowanie1 = Kp1*e1  + Ki1*dydx[4] + Kd1*pochodna1
    uchyb_poprzedni_1 = e1

    # regulator PID przegubu drugiego
    e2 = th2_ - state[2]
    pochodna2 = (e2 - uchyb_poprzedni2)/dt
    sterowanie2 = Kp2*e2 + Ki2*dydx[5] + Kd2*pochodna2
    uchyb_poprzedni_2 = e2
    
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

    dydx[1] = ((+ M22*(sterowanie1 - C12*state[3] - G11) - M12*(sterowanie2 - C21*state[1] - G21))/den)
    
    dydx[2] = state[3]

    dydx[3] = ((- M21*(sterowanie1 - C12*state[3] - G11) + M11*(sterowanie2 - C21*state[1] - G21))/den)

    dydx[4] = e1

    dydx[5] = e2
    
    return dydx


# tworzenie tablicy czasu od 0..100 próbkowanej co 0.05 sekundy
dt = 0.05
t = np.arange(0, 40, dt)

# th1 i th2 są początkowymi kątami (stopnie)
# w1 i w2 są początkowymi prędkościami kątowymi (stopnie na sekundę)
th1 = 0.0
w1 = 0.0
th2 = 0.0
w2 = 0.0

# stan początkowy
state = np.radians([th1, w1, th2, w2, e1, e2])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, state, t)

x1 = L1*sin(y[:, 0])
y1 = -L1*cos(y[:, 0])

x2 = L2*sin(y[:, 2]) + x1
y2 = -L2*cos(y[:, 2]) + y1

x1z = L1*sin(th1_r)
x2z = L2*sin(th2_r)+x1z

y1z = -L1*cos(th1_r)
y2z = -L2*cos(th2_r)+y1z


fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

zadana1, = ax.plot([], [], 'g', lw=2)
zadana1.set_data([0,x1z], [0,y1z])

zadana2, = ax.plot([], [], 'r', lw=2)
zadana2.set_data([x1z,x2z], [y1z,y2z])

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

