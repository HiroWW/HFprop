import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os 
from cycler import cycler
import colorsys
import re
import math

# load uiuc data
uiuc_data = np.loadtxt('uiuc_apc9by4.7SF_rpm5013.txt', skiprows=2)
uiuc_J = uiuc_data[:,0]
uiuc_Ct = uiuc_data[:,1]
uiuc_Cp = uiuc_data[:,2]
uiuc_eta = uiuc_data[:,3]
rpm = 5013
rps = rpm/60
R = 0.114
rho = 1.225

hfprop_Cp_smoozedtable = []
hfprop_Ct_smoozedtable = []
hfprop_eta_smoozedtable = []
qprop_Ct = []
qprop_Cp = []
qprop_eta = []

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
    subprocess.run(['python3', 'perform_analysis_refactor.py', '--airfoil', 'airfoil-qpropsame.txt'])
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
    hfprop_Cp_smoozedtable.append(Cp)
    hfprop_Ct_smoozedtable.append(Ct)
    hfprop_eta_smoozedtable.append(eta)

    os.chdir('/home/hiroaki/QPROP/myrun')
    # 実行コマンド
    command = ['../bin/qprop', 'apc9×4.7SF', 'motor-dummy', f'{V}', '5013']
    # ファイルに標準出力をリダイレクト
    with open('out.dat', 'w') as outfile:
        subprocess.run(command, stdout=outfile)
   
    with open('out.dat', 'r') as file:
        lines = file.readlines()  # Read all lines from the file

        # Check if there are at least 18 lines in the file
        if len(lines) >= 18:
            line = lines[17]  # Get the 18th line (index 17)

            # Split the line into columns
            columns = line.split()

            # Extract the specific values from the 12th, 13th, and 15th indices
            try:
                Ct = float(columns[4]) /(rho * rps**2 * (2*R)**4)
                Cq = float(columns[5]) /(rho * rps**2 * (2*R)**5)
                Cp = 2 * math.pi * Cq
                eta = float(columns[10])
            except ValueError:
                print("Error: Unable to convert data to float.")

    # Append the values to the corresponding lists
    qprop_Ct.append(Ct)
    qprop_Cp.append(Cp)
    qprop_eta.append(eta)
    os.chdir('/home/hiroaki/HFprop/ANALYSIS/')

        
os.chdir('validation')

# plot the results
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

plt.figure(dpi=300, figsize=(1.414*4, 4))
plt.plot(uiuc_J, uiuc_Cp, marker='s', linestyle='None', label='UIUC WT')
plt.plot(uiuc_J, hfprop_Cp_smoozedtable, label='inhouse')
plt.plot(uiuc_J, qprop_Cp, label='qprop')
plt.xlabel('$J$')
plt.ylabel('$C_P$')
plt.xlim(0, None)
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('Cp-J-qprop.png')
# plt.show()

plt.figure(dpi=300, figsize=(1.414*4, 4))
plt.plot(uiuc_J, uiuc_Ct, marker='s', linestyle='None', label='UIUC WT')
plt.plot(uiuc_J, hfprop_Ct_smoozedtable, label='inhouse')
plt.plot(uiuc_J, qprop_Ct, label='qprop')
plt.xlabel('$J$')
plt.ylabel('$C_T$')
plt.xlim(0, None)
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('Ct-J-qprop.png')
# plt.show()

plt.figure(dpi=300, figsize=(1.414*4, 4))
plt.plot(uiuc_J, uiuc_eta, marker='s', linestyle='None', label='UIUC WT')
plt.plot(uiuc_J, hfprop_eta_smoozedtable, label='inhouse')
plt.plot(uiuc_J, qprop_eta, label='qprop')
plt.xlabel('$J$')
plt.ylabel('$\eta$')
plt.xlim(0, None)
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('eta-J-qrop.png')
# plt.show()


plt.figure(dpi=300, figsize=(1.414*4, 4))
error = (np.array(qprop_Cp) / np.array(hfprop_Cp_smoozedtable) -1 ) * 100
plt.plot(uiuc_J, error, marker='s', label='qprop/inhouse')
plt.xlabel('$J$')
plt.ylabel('$Cp error$[\%]')
plt.xlim(0, None)
plt.ylim(0, 5)
plt.legend()
plt.tight_layout() # 余白を小さくする
plt.savefig('cperror-J-qrop.png')

