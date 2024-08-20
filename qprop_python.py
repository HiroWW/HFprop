import numpy as np
import math

# 定数
PI = math.pi
EPS = 1.0E-6
MAX_ITER = 2000

# 入力パラメータ
R = 0.1143
RPM = 5000
rho = 1.225
omega = RPM * 2 * PI / 60
B = 2
V = 8.5
CL0 = 0.1  # ここで適切な値を設定
DCLDA = 2 * PI  # ここで適切な値を設定

# r/Rの配列
n = 1000
r_R = np.linspace(0.01, 1, n)
r = r_R * R

# geometryデータの読み込みと補間
geometry = np.loadtxt('geometry.txt')
geometry_r_R = geometry[:, 0]
geometry_c_R = geometry[:, 1]
c_R = np.interp(r_R, geometry_r_R, geometry_c_R)
geometry_beta = geometry[:, 2]
beta = np.interp(r_R, geometry_r_R, geometry_beta)

# airfoilデータの読み込みと補間
alpha = np.linspace(-20, 20, n)
airfoil = np.loadtxt('airfoil.txt')
airfoil_alpha = airfoil[:, 1]
airfoil_cl = airfoil[:, 3]
airfoil_cd = airfoil[:, 2]
airfoil_cl = np.interp(alpha, airfoil_alpha, airfoil_cl)
airfoil_cd = np.interp(alpha, airfoil_alpha, airfoil_cd)

# 結果を格納する配列
psi = np.zeros(len(r_R))
gamma = np.zeros(len(r_R))
cl = np.zeros(len(r_R))

# ニュートン法の反復計算
for i in range(len(r_R)):
    Ua = V
    Ut = omega * r[i]
    U = np.sqrt(Ua**2 + Ut**2)
    PSI1 = np.arctan2(Ua, Ut)
    PSI2 = beta[i] + CL0 / DCLDA
    psi[i] = max(PSI1, PSI2)
    
    for iter in range(MAX_ITER):
        cos_psi = np.cos(psi[i])
        sin_psi = np.sin(psi[i])
        
        Wa = 0.5 * Ua + 0.5 * U * sin_psi
        Wt = 0.5 * Ut + 0.5 * U * cos_psi
        va = Wa - Ua
        vt = Ut - Wt
        aoa = beta[i] - np.arctan2(Wa, Wt)
        W = np.sqrt(Wa**2 + Wt**2)
        
        cl[i] = np.interp(aoa, alpha, airfoil_cl)
        lamda = r[i] / R * Wa / Wt
        f = B / 2 * (1 - r[i] / R) / lamda
        F = 2 / PI * np.arccos(np.clip(np.exp(-f), -1, 1))
        
        gamma[i] = vt * 4 * PI * r[i] / B * F * np.sqrt(1 + (4 * lamda * R / (PI * B * r[i]))**2)
        c = c_R[i] * R
        
        res = gamma[i] - 0.5 * c * W * cl[i]
        res_psi = -0.5 * c * (cl[i] * W * cos_psi - gamma[i] * sin_psi)  # RES_PSIに対応
        
        dpsi = -res / res_psi
        dpsi = max(-0.1, min(0.1, dpsi))
        
        # 収束判定
        if abs(dpsi) < EPS:
            break
        
        # PSIの更新
        psi[i] += dpsi
    
    else:
        print(f"Not converged at r/R = {r_R[i]} RES = {dpsi}")

# 追加の処理やプロットが必要であればここに書きます

