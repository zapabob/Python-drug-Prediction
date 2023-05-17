import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np

# CSVファイルのリスト
files = ['MPH.csv','PP.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込む
dataframes = [pd.read_csv(file) for file in files]

# 各データフレームからMPHのXLOGP3値を抽出し、リストに追加
XLOGP3_values = []
for df in dataframes:
    XLOGP3_values.extend(df[df['name'] == 'MPH']['XLOGP3'].values)

# データを正規化
XLOGP3_values = (XLOGP3_values - np.mean(XLOGP3_values)) / np.std(XLOGP3_values)

# ヒストグラムを描く
plt.hist(XLOGP3_values, bins=30, density=True)

# 正規分布のフィット
mu, std = norm.fit(XLOGP3_values)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

plt.title("Fit results: mu = %.2f,  std = %.2f" % (mu, std))
plt.show()
