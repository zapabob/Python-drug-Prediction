import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# CSVファイルのリスト
files = ['MPH.csv','PP.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込む
dataframes = [pd.read_csv(file) for file in files]

# 各データフレームから 'LOG S' と 'XLOGP3' の値を抽出し、指数関数に適用
for df in dataframes:
    df['LOG S'] = np.exp(df['LOG S'])
    df['XLOGP3'] = np.exp(df['XLOGP3'])

# プロットの準備
fig, ax = plt.subplots()

# 各データフレームに対して散布図をプロットし、回帰直線を描く
for df in dataframes:
    ax.scatter(df['XLOGP3'], df['LOG S'], alpha=0.5)
    
    # 回帰直線の計算と描画
    slope, intercept, r_value, p_value, std_err = linregress(df['XLOGP3'], df['LOG S'])
    ax.plot(df['XLOGP3'], intercept + slope*df['XLOGP3'], 'r', label=f'y={slope:.2f}x+{intercept:.2f}')

# プロットの表示
ax.set_xlabel('e^(XLOGP3)')
ax.set_ylabel('e^(LOG S)')
ax.legend()
plt.show()
