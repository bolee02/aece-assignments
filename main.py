import numpy as np
import matplotlib.pyplot as plt

dt = 3600  # 1 hour
t0 = 0
t1 = 3600*24*5  # t0 + 5 days
t2 = t1 + 3600*24*15  # t1 + 15 days
t3 = t2 + 3600*24*10  # t2 + 10 days

N_A = 6.02*10**26  # molec/kmol
rho_air = 1  # kg/m^3
NOx_init = 15*10**-12  # pptv
M_air = 29  # kg/kmol
NOx_life = 1.5*24*3600  # days
T = 298  # K
H, C, N, O = 1, 12, 14, 16
V_box = 1*10**15  # cm^3


def tokgN(C):
    kgN = C*V_box*1/N_A*N
    return kgN


X = NOx_init * rho_air * N_A * 1/(1*10**6) * 1/M_air
X_store = [X]
t = 0


def t_iterate(t_start, t_end, X_it, P_it, L_it):
    while t_start < t_end:
        X_it = (X_it + P_it*dt) * 1/(1 + L_it*dt)
        X_store.append(X_it)
        t_start += dt


L = 1/NOx_life
P = 0.12/(3600*24) * N_A/(N*V_box)  #(3.3*10**-12*np.exp(270/T))*(N+H+3*0) + (3.0*10**-12*np.exp(1500/T))*(N+4*O) + (5.1*10**-12*np.exp(210/T))*(N+3*O)
print(P)

t_iterate(t0, t1, X_store[-1], 0, L)
t_iterate(t1, t2, X_store[-1], P, L)
t_iterate(t2, t3, X_store[-1], 0, L)

time = np.arange(0, 30, 30/720)
plt.plot(time, list(map(tokgN, X_store[0: -1])))
plt.xlabel('Time (days)')
plt.ylabel('Mass of NO_x in box (kgN/box)')
plt.show()
