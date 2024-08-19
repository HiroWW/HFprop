import numpy as np
import matplotlib.pyplot as plt
import math

# general inputs
R = 0.9144
T = 869.2
RPM = 2600
omega = RPM * 2 * math.pi / 60
B = 2
V = 53.64
lbyd = 1 / 0.025

# arrays for r/R
n = 100
r_R = np.linspace(0.01, 1, n)
r = r_R * R
x = omega * r / V 
lamda = V / (omega * R)
f = B / 2 * math.sqrt(lamda**2 + 1) / lamda * (1 - r_R)
F = 2 / np.pi * np.arccos(np.exp(-f))
G = F * x**2 / (1 + x**2)

# Calculate zeta ( = v_dash / V)
I1 = 0
for i in range(n):
    I1 += 4 * G[i] * (1 - lbyd / x[i]) * r_R[i] * 1 / n

I2 = 0
for i in range(n):
    I2 += 2 * G[i] * (1 - lbyd / x[i]) / (1 + x[i]**2) * r_R[i] * 1 / n

zeta = I1 / (2 * I2) * (1 - math.sqrt(1 - 4 * I2* T / I1**2)) 
v_dash = zeta * V

# Calculate Circulation
Gamma = 2 * np.pi * r * v_dash * x / (1 + x**2) * F / B

# Plot Results
print('Advance Ratio: ', lamda)
print('I1: ', I1)
print('I2: ', I2)
# plot circulation
plt.plot(r_R, Gamma)
plt.xlabel('r/R')
plt.ylabel('Circulation')
plt.savefig('circulation.png')
plt.show()
# plot Prandlt Tip Loss Factor
plt.clf()
plt.plot(r_R, F)
plt.xlabel('r/R')
plt.ylabel('Prandlt Tip Loss Factor F')
plt.savefig('F.png')
plt.show()
# plot G
plt.clf()
plt.plot(r_R, G)
plt.xlabel('r/R')
plt.ylabel('G')
plt.savefig('G.png')

