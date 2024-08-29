import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from scipy.optimize import fsolve

PI = math.pi

class HFprop:
    def main(self):
        self.initialize()
        self.solve_psi()
        self.calculate_thrust()
        self.finalize()

    def initialize(self):
        args = self.parse_args()
        self.set_condition(args.condition)
        self.load_airfoil(args.airfoil)
        self.load_geometry(args.geometry)

    def parse_args(self):
        parser = argparse.ArgumentParser(description='BEMT propeller perform analysis tool')
        parser.add_argument('--condition', type=str, default='condition.txt', help='path to condition file')
        parser.add_argument('--geometry', type=str, default='geometry.txt', help='path to geometry file')
        parser.add_argument('--airfoil', type=str, default='airfoil.txt', help='path to airfoil file')
        return parser.parse_args()

    def set_condition(self, filename):
        self.define_r_R(25)
        self.psi = np.zeros(len(self.r_R))
        self.qpropFactorFlag = True

        conditions = self.parse_txt(filename)
        self.R = conditions['R']
        self.RPM = conditions['RPM']
        self.RHO = conditions['rho']
        self.B = conditions['B']
        self.V = conditions['V']
        self.PLOTSAVE = conditions['plotsave']
        self.OMEGA = self.RPM * 2 * PI / 60
        self.RPS = self.RPM / 60
        self.J = self.V / (self.RPS * 2 * self.R)
        self.Ua = self.V
        self.r = self.r_R * self.R

    def define_r_R(self, n):
        start = 0.15
        end = 1.0
        # get the end points of the intervals
        edges = np.linspace(start, end, n+1)
        # calculate the center of each interval
        self.r_R = (edges[:-1] + edges[1:]) / 2
    
    def parse_txt(self, filename):
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
    
    def load_airfoil(self, filename):
        self.alpha = np.linspace(-20, 20, 1000)
        airfoil = np.loadtxt(filename)
        airfoil_alpha = airfoil[:, 0]
        airfoil_cl = airfoil[:, 2]
        airfoil_cd = airfoil[:, 1]
        self.cl_alpha = np.interp(self.alpha, airfoil_alpha, airfoil_cl)
        self.cd_alpha = np.interp(self.alpha, airfoil_alpha, airfoil_cd)

    def load_geometry(self, filename):
        # load gemoetry txt and interpolate it to r/R
        geometry = np.loadtxt(filename)
        geometry_r_R = geometry[:, 0]
        geometry_c_R = geometry[:, 1]
        geometry_beta = geometry[:, 2]
        self.c_R = np.interp(self.r_R, geometry_r_R, geometry_c_R)
        self.c = self.c_R * self.R
        self.beta = np.interp(self.r_R, geometry_r_R, geometry_beta)
    
    def solve_psi(self):
        self.update()
        for i in range(len(self.r_R)):
            def equation(psi):
                self.psi[i] = psi
                self.update()
                return self.gamma[i] - 0.5 *self.W[i] * self.c[i] * self.cl[i]
            
            ans = fsolve(equation, np.arctan2(self.Ua, self.Ut[i]))
            self.psi[i] = ans
        self.update()
    
    def update(self):
        self.Ut = self.OMEGA * self.r
        self.U = np.sqrt(self.Ua**2 + self.Ut**2)
        self.Wa = 0.5 * self.Ua + 0.5 * self.U * np.sin(self.psi)
        self.Wt = 0.5 * self.Ut + 0.5 * self.U * np.cos(self.psi)
        self.phi = np.arctan2(self.Wa, self.Wt)
        self.va = self.Wa - self.Ua
        self.vt = self.Ut - self.Wt
        self.aoa = self.beta - np.degrees(np.arctan(self.Wa / self.Wt))
        self.W = np.sqrt(self.Wa**2 + self.Wt**2)
        self.cl = np.interp(self.aoa, self.alpha, self.cl_alpha)
        self.cd = np.interp(self.aoa, self.alpha, self.cd_alpha)
        self.lamda = self.r / self.R * self.Wa / self.Wt
        self.f = self.B / 2 * (1 - self.r_R) / self.lamda
        self.F = 2 / PI * np.arccos(np.clip(np.exp(-self.f), -1, 1))
        rootFactor = np.sqrt(1 + (4 * self.lamda * self.R / (PI * self.B * self.r))**2)
        self.gamma = self.vt * 4 * PI * self.r / self.B * self.F * rootFactor

    def calculate_thrust(self):
        dL = 1/2 * self.RHO * self.W**2 * self.c * self.cl * self.B
        dD = 1/2 * self.RHO * self.W**2 * self.c * self.cd * self.B
        dT = dL * np.cos(self.phi) - dD * np.sin(self.phi)
        dQ = (dL * np.sin(self.phi) + dD * np.cos(self.phi)) * self.r
        self.T = np.trapz(dT, self.r)
        self.Q = np.trapz(dQ, self.r)
        self.Ct = self.T / (self.RHO * self.RPS**2 * (2*self.R)**4)
        self.Cq = self.Q / (self.RHO * self.RPS**2 * (2*self.R)**5)
        self.Cp = self.Cq * 2 * PI
        self.eta = self.Ct * self.J / self.Cp
    
    def finalize(self):
        self.save_results()
        if (self.PLOTSAVE):
            self.plot_results()
    
    def save_results(self):
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
RADIUS = {self.R} [m]
RPM = {self.RPM}
V = {self.V} [m/s]
J = {self.J}

==== BASIC RESULTS
Thrust = {self.T} [N]
Torque = {self.Q} [Nm]
Ct = {self.Ct} [UIUC definition]
Cq = {self.Cq} [UIUC definition]
Cp = {self.Cp} [UIUC definition]
eta = {self.eta} [UIUC definition]

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
        np.savetxt('results.txt', np.array([self.r, self.cl, self.aoa, self.Wa, self.Wt, self.psi, self.Ut, self.U, self.W, self.va, self.vt, self.gamma]).T,
                header=header, fmt='%.14f')


    def plot_results(self):
        plt.plot(self.alpha, self.cl_alpha)
        plt.xlabel('alpha')
        plt.ylabel('cl')
        plt.grid()
        # plot cl
        plt.clf()
        plt.plot(self.r_R, self.cl)
        plt.xlabel('r/R')
        plt.ylabel('Cl')
        plt.grid()
        plt.show()
        # plot aoa
        plt.clf()
        plt.plot(self.r_R, self.aoa)
        plt.xlabel('r/R')
        plt.ylabel('aoa')
        plt.grid()
        # plt.show()

# ============================================== RUN THE COE ============================================== 
hfprop = HFprop()
hfprop.main()
# =========================================================================================================