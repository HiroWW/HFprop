import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os 

# load uiuc data
uiuc_data = np.loadtxt('uiuc_apc9by4.7SF_rpm5013.txt', skiprows=2)
uiuc_J = uiuc_data[:,0]
uiuc_Ct = uiuc_data[:,1]
uiuc_Cp = uiuc_data[:,2]
uiuc_eta = uiuc_data[:,3]
rpm = 5013
rps = rpm/60
R = 0.1143

hfprop_Cp = []
hfprop_Ct = []
hfprop_eta = []

os.chdir('..')

for J in uiuc_J:
    V = J * rps * 2 * R
    # refresh and write condition.txt fil
    with open('condition.txt', 'w') as f:
        f.write('# condition for CALCULATION\n')
        f.write('R = 0.1143\n')
        f.write('RPM = {}\n'.format(rpm))
        f.write('rho = 1.225\n')
        f.write('B = 2\n')
        f.write('V = {}\n'.format(V))
        f.write('plotsave = False\n')
    # run the calculation
    subprocess.run(['python3', 'perform_analysis.py'])
    # load the results and append to the list
    with open('results.txt', 'r') as file:
        for line in file:
            # Skip lines that don't contain the relevant data
            if not line.startswith('#'):
                continue

            # Check and extract the values for Ct, Cp, and eta
            if 'Ct =' in line:
                Ct = float(line.split('=')[1].strip().split()[0])
            elif 'Cp =' in line:
                Cp = float(line.split('=')[1].strip().split()[0])
            elif 'eta =' in line:
                eta = float(line.split('=')[1].strip().split()[0])
    hfprop_Cp.append(Cp)
    hfprop_Ct.append(Ct)
    hfprop_eta.append(eta)
        
os.chdir('validation')

# plot the results
plt.plot(uiuc_J, uiuc_Cp, label='UIUC')
plt.plot(uiuc_J, hfprop_Cp, 'o', label='HFPROP')
plt.xlabel('J')
plt.ylabel('Cp')
plt.legend()
plt.savefig('Cp-J.png')
plt.show()
plt.clf()
plt.plot(uiuc_J, uiuc_Ct, label='UIUC')
plt.plot(uiuc_J, hfprop_Ct, 'o', label='HFPROP')
plt.xlabel('J')
plt.ylabel('Ct')
plt.legend()
plt.savefig('Ct-J.png')
plt.show()
plt.clf()
plt.plot(uiuc_J, uiuc_eta, label='UIUC')
plt.plot(uiuc_J, hfprop_eta, 'o', label='HFPROP')
plt.xlabel('J')
plt.ylabel('eta')
plt.legend()
plt.savefig('eta-J.png')
plt.show()
