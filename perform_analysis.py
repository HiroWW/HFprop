import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve

# general inputs
# for now, use the same values as in Kawasaki et.al, 2024 to validate the code
R = 0.1143
RPM = 5000
rho = 1.225
omega = RPM * 2 * math.pi / 60
B = 2
V = 8.5
# dbyl = 0.025
# lbyd = 1 / dbyl
# Tc = 2 * T / (rho * V**2 * R**2 * math.pi)

# arrays for r/R
n = 1000
r_R = np.linspace(0.01, 1, n)
r = r_R * R
x = omega * r / V 
lamda = V / (omega * R)
f = B / 2 * math.sqrt(lamda**2 + 1) / lamda * (1 - r_R)
F = 2 / np.pi * np.arccos(np.exp(-f))
G = F * x**2 / (1 + x**2)

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

# def equation(v_dash):
#     aoa = beta - np.arctan(x * R / r * (1 + 1/2 * v_dash / V))
#     cl = np.interp(aoa, alpha, airfoil_cl)   
#     return V * c_R * cl / v_dash - 4 * math.pi / B * lamda / np.sqrt(1 + x**2) * G

# v_dash = fsolve(equation, 1)

# v_dash = np.zeros_like(r_R)  # v_dashを初期化

# for i in range(len(r_R)):

def equation(v_dash_i):
    i = n - 2
    aoa = beta[i] - np.arctan(x[i] * R / r[i] * (1 + 1/2 * v_dash_i / V))
    cl = np.interp(aoa, alpha, airfoil_cl)
    return V * c_R[i] * cl / v_dash_i - 4 * math.pi / B * lamda / np.sqrt(1 + x[i]**2) * G[i]

v_dash = fsolve(equation, 1)

# Calculate Circulation
Gamma = 2 * np.pi * r * v_dash * x / (1 + x**2) * F / B

# Plot Results
print('Advance Ratio: ', lamda)
# plot circulation
plt.plot(r_R, Gamma)
plt.xlabel('r/R')
plt.ylabel('Circulation')
plt.savefig('circulation.png')
# plot v_dash
plt.clf()
# plt.plot(r_R, v_dash)
# # set y-axis limit
# plt.ylim(0, 100.5)
# plt.xlabel('r/R')
# plt.ylabel('v_dash')
# plt.savefig('v_dash.png')
