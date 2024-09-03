import numpy as np
import matplotlib.pyplot as plt
import math
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os 
from cycler import cycler
import colorsys
import re
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
# =========================================
# ========= Plotting Congifuration ========
# =========================================
# アスペクト比を大和比に設定
# 参考：https://qiita.com/ShotaDeguchi/items/c24a89272e651c3de58c
# Font settings
# 参考：https://qiita.com/yuki_2020/items/7f96c79614c2576a37d3
# 参考；https://qiita.com/nabenabe0928/items/a026555917525e2bd761(timesのinstall)
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'Times New Roman' 
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams["font.size"] = 15 
# 色をはっきりしたものに変更
# 参考：https://maskit.hatenablog.com/entry/2022/05/29/124006
clist = []
for i in range(8):
    c = colorsys.rgb_to_hsv( *plt.get_cmap('Set1')(i)[:3] )
    c = colorsys.hsv_to_rgb( c[0] , 1.-(1.-c[1])*.2 , c[2]*.85 )
    clist.append( c )
plt.rcParams['axes.prop_cycle'] = cycler( color=clist )
# 軸、グリッド、目盛りの設定
# 参考：https://qiita.com/yuki_2020/items/7f96c79614c2576a37d3
plt.rcParams['xtick.direction'] = 'in' #x軸の目盛りの向き
plt.rcParams['ytick.direction'] = 'in' #y軸の目盛りの向き
plt.rcParams['axes.grid'] = True # グリッドの作成
plt.rcParams['grid.linestyle']='-' #グリッドの線種
plt.rcParams["xtick.minor.visible"] = True  #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True  #y軸補助目盛りの追加
plt.rcParams['xtick.top'] = True  #x軸の上部目盛り
plt.rcParams['ytick.right'] = True  #y軸の右部目盛り
plt.rcParams['axes.linewidth'] = 1.0 # 軸の太さ
# markerの設定
# 参考：https://www.useful-python.com/matplotlib-thesis-format/
plt.rcParams['lines.markerfacecolor'] = 'white' #中を白塗り
plt.rcParams['lines.markersize'] = 6 #マーカーサイズ
# 凡例の設定
plt.rcParams["legend.fancybox"] = False  # 丸角OFF
plt.rcParams["legend.framealpha"] = 1  # 透明度の指定、0で塗りつぶしなし
plt.rcParams["legend.edgecolor"] = 'black'  # edgeの色を変更
# plt.rcParams["legend.markerscale"] = 3 #markerサイズの倍率
# cl
plt.figure(dpi=300, figsize=(1.414*4, 4))
plt.plot(alpha, cl, label='utcart')
plt.plot(alpha, cl_q, label='qprop')
plt.xlabel('$\\alpha$')
plt.ylabel('$C_l$')
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('cl-alpha-compare.png')

# plt.show()
# cd
plt.figure(dpi=300, figsize=(1.414*4, 4))
plt.plot(alpha, cd, label='utcart')
plt.plot(alpha, cd_q, label='qprop')
plt.xlabel('$\\alpha$')
plt.ylabel('$C_d$')
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('cd-alpha-compare.png')
# plt.show()