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
# dbyl = 0.025
# lbyd = 1 / dbyl
# Tc = 2 * T / (rho * V**2 * R**2 * math.pi)

# arrays for r/R
n = 1000
r_R = np.linspace(0.01, 1, n)
r = r_R * R
lamda = V / (omega * R)
f = B / 2 * math.sqrt(lamda**2 + 1) / lamda * (1 - r_R)
F = 2 / np.pi * np.arccos(np.exp(-f))

# x = omega * r / V 
# G = F * x**2 / (1 + x**2)

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

# r_Rと同じ要素数の二次元配列solutionを作成
solution = np.zeros((len(r_R), 2))

# r_Rと同じ要素数の二次元配列solutionを作成
solution = np.zeros((len(r_R), 2))
# solu = np.zeros_like(r_R)  # v_dashを初期化

for i in range(len(r_R)):

    def equations(vars):
        phi_b, phi_s = vars
        w =  (math.sin(phi_b) - V) / math.cos(phi_b)
        aoa = beta[i] - phi_b
        cl = np.interp(aoa, alpha, airfoil_cl)
        W = (w*math.cos(phi_b) + V) / math.cos(phi_b)
        dl = 1/2 * rho * W**2 * c_R[i] * cl
        eq1 = (math.sin(phi_b) - V) / math.cos(phi_b) - (math.sin(phi_s) - V) / math.cos(phi_s)
        eq2 = rho * (omega*r_R[i]*R - w*math.sin(phi_b)) / math.cos(phi_b) * 2 * math.pi / B * r_R[i] * R * 0.5*w * math.sin(phi_s) * F[i] - dl
        return [eq1, eq2]

    solution[i] = fsolve(equations, (1, 1))

# convert to degree and print solution
phi_b = solution[:, 0] * 180 / math.pi
phi_s = solution[:, 1] * 180 / math.pi

print('phi_b: ', phi_b[0])
print('phi_s: ', phi_s[0])
#plot phi_b and phi_s
plt.plot(r_R, phi_b)
plt.plot(r_R, phi_s)
plt.xlabel('r/R')
plt.ylabel('phi_b and phi_s')
plt.savefig('phi_b_s.png')


# # Calculate Circulation
# Gamma = 2 * np.pi * r * v_dash * x / (1 + x**2) * F / B

# # Calculate induced velocity component tangential to the rotor plane
# print('v_dash: ', v_dash)
# v_i = v_dash * x**2 / (1 + x**2)

# # Plot Results
# print('Advance Ratio: ', lamda)
# # plot circulation
# plt.plot(r_R, Gamma)
# plt.xlabel('r/R')
# plt.ylabel('Circulation')
# plt.savefig('circulation.png')
# # plot v_dash
# plt.clf()
# plt.plot(r_R, v_i/340/2*F)
# # set y-axis limit
# plt.ylim(0, 0.2)
# plt.xlabel('r/R')
# plt.ylabel('v_i')
# plt.savefig('v_i.png')
