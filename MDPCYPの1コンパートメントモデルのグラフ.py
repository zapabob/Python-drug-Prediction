import numpy as np
import matplotlib.pyplot as plt
import os

def one_compartment_model(dose_mg, BA, weight_kg, Vd_L, half_life_h, t_end=24, num_points=1000):
    k_el = np.log(2) / half_life_h  # 排泄定数（h^-1）

    # 実際に吸収される投与量（mg）
    dose_mg_absorbed = dose_mg * BA

    # 初期濃度（mg/L）
    C0_mg_per_L = dose_mg_absorbed / Vd_L

    # 時間範囲を設定
    t = np.linspace(0, t_end, num_points)  # 時間（h）

    # 1コンパートメントモデルの濃度-時間曲線を計算
    C_mg_per_L = C0_mg_per_L * np.exp(-k_el * t)
    
    return t, C_mg_per_L

def plot_and_save(t, C_mg_per_L):
    plt.figure(figsize=(10, 6))
    plt.plot(t, C_mg_per_L)
    plt.xlabel('Time (h)')
    plt.ylabel('Drug concentration (mg/L)')
    plt.title('One-compartment model: concentration vs time')
    plt.grid(True)
    # デスクトップに保存する
    plt.savefig(os.path.expanduser('~/Desktop/one_compartment_model.png'))
    plt.show()

if __name__ == "__main__":
    dose_mg = 15
    BA = 0.98  # バイオアベイラビリティ（98%）
    weight_kg = 60
    Vd_L = 36
    half_life_h = 4.5
    
    t, C_mg_per_L = one_compartment_model(dose_mg, BA, weight_kg, Vd_L, half_life_h)
    plot_and_save(t, C_mg_per_L)

