import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz  # AUCを求めるための関数

def one_compartment_model(dose_mg, Vd_L, half_life_h, MW, t_end=24, num_points=1000):
    k_el = np.log(2) / half_life_h  # 排泄定数（h^-1）

    # 投与量をモル単位に変換
    dose_mol = dose_mg / MW * 1e-3  # モル単位への変換（mol）

    # 初期濃度（モル/L）
    C0_mol_per_L = dose_mol / Vd_L

    # 時間範囲を設定
    t = np.linspace(0, t_end, num_points)  # 時間（h）

    # 1コンパートメントモデルの濃度-時間曲線を計算
    C_mol_per_L = C0_mol_per_L * np.exp(-k_el * t)
    
    return t, C_mol_per_L

def plot_and_save(t, C_mol_per_L, filename='one_compartment_model.png'):
    plt.figure(figsize=(10, 6))
    plt.plot(t, C_mol_per_L)
    plt.xlabel('Time (h)')
    plt.ylabel('Drug concentration (mol/L)')
    plt.title('One-compartment model: concentration vs time')
    plt.grid(True)
    plt.savefig(filename)
    plt.show()

if __name__ == "__main__":
    dose_mg = 15
    Vd_L = 36
    half_life_h = 4.5
    MW = 315.41
    
    t, C_mol_per_L = one_compartment_model(dose_mg, Vd_L, half_life_h, MW)
    plot_and_save(t, C_mol_per_L, filename='C:/Users/YourUsername/Desktop/one_compartment_model.png')
