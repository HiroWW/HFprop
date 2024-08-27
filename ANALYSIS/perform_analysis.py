import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from scipy.optimize import fsolve
from scipy.optimize import minimize

# parse arguments
def parse_txt(filename):
    txt_data = {}
    with open(filename, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith('#') or not line.strip():
                continue
            # Split the line into key and value and store them in the dictionary
            key, value = line.split('=')
            key = key.strip()
            value = value.strip()
            # Handle boolean values
            if value.lower() == 'true':
                txt_data[key] = True
            elif value.lower() == 'false':
                txt_data[key] = False
            else:
                try:
                    # Attempt to convert to a numeric value
                    txt_data[key] = float(value)
                except ValueError:
                    # If not numeric, store as a string
                    txt_data[key] = value
    return txt_data
parser = argparse.ArgumentParser(description='BEMT propeller perform analysis tool')
parser.add_argument('--condition', type=str, default='condition.txt', help='path to condition file')
parser.add_argument('--geometry', type=str, default='geometry.txt', help='path to geometry file')
parser.add_argument('--airfoil', type=str, default='airfoil.txt', help='path to airfoil file')
args = parser.parse_args()

# general inputs
conditions = parse_txt(args.condition)
R = conditions['R']
RPM = conditions['RPM']
rho = conditions['rho']
B = conditions['B']
V = conditions['V']
plotsave = conditions['plotsave']

omega = RPM * 2 * math.pi / 60
rps = RPM / 60
J = V / (rps * 2 * R)
CL0 = math.radians(-6)  # example
DCLDA = 2 * math.pi  # example

# arrays for r/R
n = 25
start = 0.15
end = 1.0
# get the end points of the intervals
edges = np.linspace(start, end, n+1)
# calculate the center of each interval
r_R = (edges[:-1] + edges[1:]) / 2
r = r_R * R

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
airfoil_alpha = airfoil[:, 0]
airfoil_cl = airfoil[:, 2]
airfoil_cd = airfoil[:, 1]
airfoil_cl = np.interp(alpha, airfoil_alpha, airfoil_cl)
airfoil_cd = np.interp(alpha, airfoil_alpha, airfoil_cd)

plt.plot(alpha, airfoil_cl)
plt.xlabel('alpha')
plt.ylabel('cl')
plt.savefig('./debug/cl-alpha.png')
plt.clf()
plt.plot(alpha, airfoil_cd)
plt.xlabel('alpha')
plt.ylabel('cd')
plt.savefig('./debug/cd-alpha.png')

plt.clf()

# solve for psi by newton method
psi = np.zeros(len(r_R))
qpropFactorFlag = True
for i in range(len(r_R)):
    def equation(psi):
        r = r_R[i] * R
        Ua = V
        Ut = omega * r
        U = math.sqrt(Ua**2 + Ut**2)
        Wa = 1/2 * Ua + 1/2 * U * np.sin(psi)
        Wt = 1/2 * Ut + 1/2 * U * np.cos(psi)
        va = Wa - Ua
        vt = Ut - Wt
        aoa = beta[i] - np.degrees(np.arctan(Wa / Wt))
        W = np.sqrt(Wa**2 + Wt**2)
        cl = np.interp(aoa, alpha, airfoil_cl) * 1 / np.sqrt(1-(W/340)**2)
        # --- qprop modify ---
        if (qpropFactorFlag):
            lamda = r / R * Wa / Wt
            f = B / 2 * (1 - r/R) / lamda 
            f = np.minimum(f, 20)
            F = 2 / math.pi * np.arccos(np.clip(np.exp(-f), -1, 1))
            gamma = vt * 4 * math.pi * r / B * F * np.sqrt(1 + (4 * lamda * R / (math.pi * B * r ))**2)
        # --- prandtl original ---
        else:
            lamda = V / R / omega 
            f = B / 2 * (1 - r/R) / lamda *np.sqrt(lamda**2 + 1)
            F = 2 / math.pi * np.arccos(np.clip(np.exp(-f), -1, 1))
            gamma = vt * 4 * math.pi * r / B * F
        c = c_R[i] * R
        return gamma - 1/2 * W * c * cl
    # initial guess for psi (from QPROP, Drela, 2007)
    Ua = V
    Ut = omega * r[i]
    initial_guess = np.maximum(np.arctan2(Ua, Ut), math.radians(beta[i])+ CL0 / DCLDA)
    psi[i] = fsolve(equation, initial_guess)[0]

# plot psi
plt.plot(r_R, psi)
plt.xlabel('r/R')
plt.ylabel('psi')
plt.savefig('./debug/psi.png')

# post process
r = r_R * R
c = c_R * R
Ua = V
Ut = omega * r
U = np.sqrt(Ua**2 + Ut**2)
Wa = 1/2 * Ua + 1/2 * U * np.sin(psi)
Wt = 1/2 * Ut + 1/2 * U * np.cos(psi)
W = np.sqrt(Wa**2 + Wt**2)
va = Wa - Ua
vt = Ut - Wt
phi = np.arctan(Wa / Wt)
aoa = beta - np.degrees(phi)
cl = np.interp(aoa, alpha, airfoil_cl)
cd = np.interp(aoa, alpha, airfoil_cd)
dL = 1/2 * rho * W**2 * c * cl * B
dD = 1/2 * rho * W**2 * c * cd * B
dT = dL * np.cos(phi) - dD * np.sin(phi)
dQ = (dL * np.sin(phi) + dD * np.cos(phi)) * r
T = np.trapz(dT, r)
Q = np.trapz(dQ, r)
Ct = T / (rho * rps**2 * (2*R)**4)
Cq = Q / (rho * rps**2 * (2*R)**5)
Cp = Cq * 2 * math.pi
if qpropFactorFlag:
    lamda = r / R * Wa / Wt
    f = B / 2 * (1 - r/R) / lamda
    F = 2 / math.pi * np.arccos(np.clip(np.exp(-f), -1, 1))
    gamma = vt * 4 * math.pi * r / B * F * np.sqrt(1 + (4 * lamda * R / (math.pi * B * r ))**2)
else:
    lamda = V / R / omega 
    f = B / 2 * (1 - r/R) / lamda *np.sqrt(lamda**2 + 1)
    F = 2 / math.pi * np.arccos(np.clip(np.exp(-f), -1, 1))
    gamma = vt * 4 * math.pi * r / B * F

# Define the header with explanations
header = f'''HFPROP RESULT
██╗  ██╗███████╗██████╗ ██████╗  ██████╗ ██████╗ 
██║  ██║██╔════╝██╔══██╗██╔══██╗██╔═══██╗██╔══██╗
███████║█████╗  ██████╔╝██████╔╝██║   ██║██████╔╝
██╔══██║██╔══╝  ██╔═══╝ ██╔══██╗██║   ██║██╔═══╝ 
██║  ██║██║     ██║     ██║  ██║╚██████╔╝██║     
╚═╝  ╚═╝╚═╝     ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝

Hiroaki Fujiwara, 2024

==== CALCULATION CONDITIONS
RADIUS = {R} [m]
RPM = {RPM}
V = {V} [m/s]
J = {J}

==== BASIC RESULTS
Thrust = {T} [N]
Torque = {Q} [Nm]
Ct = {Ct} [UIUC definition]
Cq = {Cq} [UIUC definition]
Cp = {Cp} [UIUC definition]
eta = {Ct * J / Cp} [UIUC definition]

==== BLADE ELEMENT RESULTS
r: radial position [m]
cl: lift coefficient [-]
aoa: angle of attack [deg]
Wa: axial velocity component [m/s]
Wt: tangential velocity component [m/s]
psi: azimuthal position [rad]
Ut: tangential velocity [m/s]
U: total velocity [m/s]
W: resultant velocity [m/s]
va: induced axial velocity [m/s]
vt: induced tangential velocity [m/s]
gamma: circulation [-]
r              cl                aoa                Wa              Wt                psi               Ut                U                 W               va                vt                  gamma
'''

# Save the results to 'results.txt' with the specified format
np.savetxt('results.txt', np.array([r, cl, aoa, Wa, Wt, psi, Ut, U, W, va, vt, gamma]).T,
           header=header, fmt='%.14f')

# load qprop data for validation
qprop = np.loadtxt('../secret/qprop_validation_advance.txt')
qprop_r_R = qprop[:, 0] / R
qprop_c_R = qprop[:, 1] / R
qprop_beta = qprop[:, 2]
qprop_Wa  = qprop[:, 9]
qprop_cl  = qprop[:, 3]
qprop_psi = np.arcsin((qprop_Wa - 0.5*V) / (0.5*U))
qprop_Wt = 0.5*Ut + 0.5*U*np.cos(qprop_psi)
qprop_W = np.sqrt(qprop_Wa**2 + qprop_Wt**2) # Mach数を使わないで計算。
qprop_lamda = r / R * qprop_Wa / qprop_Wt
qprop_f = B / 2 * (1 - r/R) / qprop_lamda
qprop_F = 2 / math.pi * np.arccos(np.clip(np.exp(-qprop_f), -1, 1))
qprop_vt = omega * r - qprop_Wt
qprop_r = qprop_r_R * R
qprop_Ut = omega * qprop_r
qprop_U = np.sqrt(V**2 + qprop_Ut**2)
qprop_gamma = qprop_vt * 4 * math.pi * r / B * qprop_F * np.sqrt(1 + (4 * qprop_lamda * R / (math.pi * B * r ))**2)

if (plotsave):
    # plot results
    # plot F
    plt.clf()
    plt.plot(r_R, F, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_F, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('F')
    plt.savefig('./debug/F.png')
    # plot vt
    plt.clf()
    plt.plot(r_R, vt, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_vt, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('vt')
    plt.savefig('./debug/vt.png')
    # plot Wt
    plt.clf()
    plt.plot(r_R, qprop_Wt, marker='o',label='in-house')
    plt.plot(qprop_r_R, 0.5*qprop_Ut+0.5*qprop_U*np.cos(qprop_psi), marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('Wt')
    plt.savefig('./debug/Wt.png')
    # plot c_R
    plt.clf()
    plt.plot(r_R, c_R, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_c_R, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('c_R')
    plt.savefig('./debug/c_R.png')
    # plot beta
    plt.clf()
    plt.plot(r_R, beta, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_beta, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('beta')
    plt.savefig('./debug/beta.png')
    # plot W
    plt.clf()
    plt.plot(r_R, W, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_W, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('W')
    plt.savefig('./debug/W.png')
    # plot psi
    plt.clf()
    plt.plot(r_R, psi, marker='o',label='in-house')
    plt.plot(qprop_r_R, qprop_psi, marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('psi')
    plt.savefig('./debug/psi.png')
    # plot gamma
    plt.clf()
    plt.plot(r_R, gamma, marker='o',label='in-house(wake)')
    plt.plot(r_R, 0.5*W*c_R*R*cl, marker='o',label='in-house(on-blade)')
    plt.plot(qprop_r_R, 0.5*qprop_W*qprop_c_R*R*qprop_cl, marker='o',label='qprop(on-blade)')
    plt.plot(qprop_r_R, qprop_gamma, marker='o',label='qprop(wake)')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('gamma')
    plt.savefig('./debug/gamma.png')
    # plot va
    plt.clf()
    plt.plot(r_R, -va/340*F)
    plt.xlabel('r/R')
    plt.ylabel('va')
    plt.savefig('./debug/va.png')
    # plot Wa
    plt.clf()
    plt.plot(r_R, Wa, marker='o', label='in-house')
    plt.plot(qprop_r_R, qprop_Wa,  marker='o',label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('Wa')
    plt.savefig('./debug/Wa.png')
    #plot cl
    plt.clf()
    plt.plot(r_R, cl, marker='o', label='in-house')
    plt.plot(qprop_r_R, qprop_cl, marker='o', label='qprop')
    plt.legend()
    plt.xlabel('r/R')
    plt.ylabel('cl')
    plt.savefig('./debug/cl-r.png')
