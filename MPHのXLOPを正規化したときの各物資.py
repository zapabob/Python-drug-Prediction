import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np

# CSVファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込む
dataframes = [pd.read_csv(file) for file in files]

# 各データフレームからMPHのXLOGP値を抽出し、リストに追加
xlogp_values = []
for df in dataframes:
    xlogp_values.extend(df[df['name'] == 'MPH']['XLOGP'].values)

# データを正規化
xlogp_values = (xlogp_values - np.mean(xlogp_values)) / np.std(xlogp_values)

# ヒストグラムを描く
plt.hist(xlogp_values, bins=30, density=True)

# 正規分布のフィット
mu, std = norm.fit(xlogp_values)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

plt.title("Fit results: mu = %.2f,  std = %.2f" % (mu, std))
plt.show()
