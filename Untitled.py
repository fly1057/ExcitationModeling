import numpy as np
import matplotlib.pyplot as plt


def rk4(y, f, dt, t):
    k1 = f(t, y)
    k2 = f(t + dt / 2, y + dt / 2 * k1)
    k3 = f(t + dt / 2, y + dt / 2 * k2)
    k4 = f(t + dt, y + dt * k3)
    ynext = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return ynext


def f(t, y):
    #y=t^3-t^2+3
    #y'=3t^2-2t
    return 3 * t * t - 2 * t


dt = 0.01
t0 = 0.0
y0 = 3.0
t = t0
y = y0
tseq = []
yseq = []
tseq_original = np.arange(0,1,0.01)
yseq_original = tseq_original**3-tseq_original**2+3
while t < 1:
    tseq.append(t)
    yseq.append(y)
    y = rk4(y, f, dt, t)
    t = t + dt

plt.plot(tseq, yseq,'+-')
plt.plot(tseq_original,yseq_original,'-')
plt.show()
