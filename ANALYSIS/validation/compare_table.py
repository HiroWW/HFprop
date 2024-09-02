import numpy as np
import matplotlib.pyplot as plt
import math

# utcart
filename = "../airfoil.txt"
airfoil = np.loadtxt(filename)
airfoil_alpha = airfoil[:, 0]
airfoil_cl = airfoil[:, 2]
airfoil_cd = airfoil[:, 1]

alpha = np.linspace(-20, 20, 500)
cl = np.interp(alpha, airfoil_alpha, airfoil_cl)
cd = np.interp(alpha, airfoil_alpha, airfoil_cd)

# qprop way
# parameters for generating the airfoil polars
cl0 = 0.58453281
cla = 5.792954320080275
dclda = cla
cl_min = -0.25232311
cl_max = 1.4911173
cd0 = 0.020439062369484066
cd2u = 0.08436756949113908 
cd2l = 0.06754529977189071
clcd0 = 0.8051186313378577
ReynoldsEff = 70262
REexp = -0.5

cl_q = []
cd_q = []
for aoa in alpha:
    a = math.radians(aoa)
    cl_qprop = dclda*a + cl0
    # if stall, cl dicreases
    stall = False
    if (cl_qprop < cl_min):
        stall = True
        cl_qprop = cl_min * math.cos(a) - cl0/cla
    elif (cl_qprop > cl_max):
        stall = True
        cl_qprop = cl_max * math.cos(a) - cl0/cla
    
    # set cd slope based on cl
    if (cl_qprop >= clcd0):
        cd2 = cd2u
    if (cl_qprop < clcd0):
        cd2 = cd2l
    cd_qprop = cd0 + cd2*(cl_qprop-clcd0)**2
    
    # if stall, cd increases
    if (stall):
        acd0 = (clcd0-cl0)/dclda
        dcd = 2.0 * math.sin(a-acd0)**2
        cd_qprop = cd_qprop + dcd
    
    # append to array
    cl_q.append(cl_qprop)
    cd_q.append(cd_qprop)

# plot results
# cl
plt.plot(alpha, cl, label='utcart')
plt.plot(alpha, cl_q, label='qprop')
plt.legend()
plt.savefig('cl-alpha-compare.png')
# plt.show()
# cd
plt.clf()
plt.plot(alpha, cd, label='utcart')
plt.plot(alpha, cd_q, label='qprop')
plt.legend()
plt.savefig('cd-alpha-compare.png')
# plt.show()