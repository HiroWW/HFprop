import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
from scipy.optimize import minimize
# general inputs
# for now, use the same values as in Kawasaki et.al, 2024 to validate the code
R = 0.1143
RPM = 5000
rho = 1.225
omega = RPM * 2 * math.pi / 60
B = 2
V = 8.5

# arrays for r/R
n = 1000
r_R = np.linspace(0.01, 1, n)
r = r_R * R
lamda = V / (omega * R)
f = B / 2 * math.sqrt(lamda**2 + 1) / lamda * (1 - r_R)
F = 2 / np.pi * np.arccos(np.exp(-f))

# load gemoetry txt and interpolate it to r/R
geometry = np.loadtxt('geometry.txt')
geometry_r_R = geometry[:, 0]
geometry_c_R = geometry[:, 1]
c_R = np.interp(r_R, geometry_r_R, geometry_c_R)
geometry_beta = geometry[:, 2]
beta = np.interp(r_R, geometry_r_R, geometry_beta)

# load airfoil aerodynamic data and interpolate it to alpha
alpha = np.linspace(-20, 20, n)
airfoil = np.loadtxt('airfoil.txt')
airfoil_alpha = airfoil[:, 1]
airfoil_cl = airfoil[:, 3]
airfoil_cd = airfoil[:, 2]
airfoil_cl = np.interp(alpha, airfoil_alpha, airfoil_cl)
airfoil_cd = np.interp(alpha, airfoil_alpha, airfoil_cd)


psis = np.zeros(len(r_R))

for i in range(len(r_R)):
    def equation(psi):
        r = r_R[i] * R
        Ua = V
        Ut = omega * r
        U = math.sqrt(Ua**2 + Ut**2)
        Wa = 1/2 * Ua + 1/2 * U * np.sin(psi)
        Wt = 1/2 * Ut + 1/2 * U * np.cos(psi)
        W = np.sqrt(Wa**2 + Wt**2)
        va = Wa - Ua
        vt = Ut - Wt
        aoa = beta[i] - np.arctan(Wa / Wt)
        cl = np.interp(aoa, alpha, airfoil_cl)
        gamma = vt * 4 * math.pi * r / B * F[i]
        c = c_R[i] * R
        return gamma - 1/2 * W * c * cl

    psis[i] = fsolve(equation, 1)[0]

# plot psi
plt.plot(r_R, psis)
plt.xlabel('r/R')
plt.ylabel('psi')
plt.savefig('./debug/psi.png')


r = r_R * R
Ua = V
Ut = omega * r
U = np.sqrt(Ua**2 + Ut**2)
Wa = 1/2 * Ua + 1/2 * U * np.sin(psis)
Wt = 1/2 * Ut + 1/2 * U * np.cos(psis)
W = np.sqrt(Wa**2 + Wt**2)
va = Wa - Ua
vt = Ut - Wt
aoa = beta - np.arctan(Wa / Wt)
cl = np.interp(aoa, alpha, airfoil_cl)
gamma = vt * 4 * math.pi * r / B * F

# plot gamma
plt.clf()
plt.plot(r_R, -va/340)
plt.xlabel('r/R')
plt.ylabel('va')
plt.savefig('./debug/va.png')
