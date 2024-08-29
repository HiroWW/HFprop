import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from scipy.optimize import fsolve

PI = math.pi

class HFprop:    
    def initialize(self):

        self.define_r_R(25)

        parser = argparse.ArgumentParser(description='BEMT propeller perform analysis tool')
        parser.add_argument('--condition', type=str, default='condition.txt', help='path to condition file')
        parser.add_argument('--geometry', type=str, default='geometry.txt', help='path to geometry file')
        parser.add_argument('--airfoil', type=str, default='airfoil.txt', help='path to airfoil file')
        args = parser.parse_args()
        
        self.load_condition(args.condition)
        
        self.OMEGA = self.RPM * 2 * PI / 60
        self.RPS = self.RPM / 60
        self.J = self.V / (self.RPS * 2 * self.R)
        self.Ua = self.V
        self.r = self.r_R * self.R

        self.load_airfoil(args.airfoil)

        self.load_geometry(args.geometry)

        self.psi = np.zeros(len(self.r_R))
        self.qpropFactorFlag = True

    def load_condition(self, filename):
        conditions = self.parse_txt(filename)
        self.R = conditions['R']
        self.RPM = conditions['RPM']
        self.RHO = conditions['rho']
        self.B = conditions['B']
        self.V = conditions['V']
        self.PLOTSAVE = conditions['plotsave']
    
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
        self.cl_alpha = np.interp(self.alpha, airfoil_alpha, airfoil_cl)

    
    def load_geometry(self, filename):
        # load gemoetry txt and interpolate it to r/R
        geometry = np.loadtxt(filename)
        geometry_r_R = geometry[:, 0]
        geometry_c_R = geometry[:, 1]
        geometry_beta = geometry[:, 2]
        self.c_R = np.interp(self.r_R, geometry_r_R, geometry_c_R)
        self.c = self.c_R * self.R
        self.beta = np.interp(self.r_R, geometry_r_R, geometry_beta)
    
    def define_r_R(self, n):
        start = 0.15
        end = 1.0
        # get the end points of the intervals
        edges = np.linspace(start, end, n+1)
        # calculate the center of each interval
        self.r_R = (edges[:-1] + edges[1:]) / 2
        
    def update(self):
        self.Ut = self.OMEGA * self.r
        self.U = np.sqrt(self.Ua**2 + self.Ut**2)
        self.Wa = 0.5 * self.Ua + 0.5 * self.U * np.sin(self.psi)
        self.Wt = 0.5 * self.Ut + 0.5 * self.U * np.cos(self.psi)
        self.va = self.Wa - self.Ua
        self.vt = self.Ut - self.Wt
        self.aoa = self.beta - np.degrees(np.arctan(self.Wa / self.Wt))
        self.W = np.sqrt(self.Wa**2 + self.Wt**2)
        self.cl = np.interp(self.aoa, self.alpha, self.cl_alpha)
        # self.lamda = self.r / self.R * self.Wa / self.Wt
        # self.f = self.B / 2 * (1 - self.r_R) / self.lamda
        # self.F = 2 / PI * np.arccos(np.clip(np.exp(-self.f), -1, 1))
        # rootFactor = np.sqrt(1 + (4 * self.lamda * self.R / (PI * self.B * self.r))**2)
        self.gamma = self.vt * 4 * PI * self.r / self.B #* self.F * rootFactor

    def solve_psi(self):
        self.update()
        for i in range(len(self.r_R)):
            def equation(psi):
                self.psi[i] = psi
                self.update()
                return self.gamma[i] - 0.5 *self.W[i] * self.c[i] * self.cl[i]
            
            ans = fsolve(equation, np.arctan2(self.Ua, self.Ut[i]))
            self.psi[i] = ans
            print("ans:", self.psi[i])
    
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

    
    def main(self):
        self.initialize()
        self.solve_psi()
        self.update()
        self.plot_results()


# ============================================== RUN THE COE ============================================== 
hfprop = HFprop()
hfprop.main()
# =========================================================================================================