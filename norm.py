import numpy as np
import pandas as pd
import matplotlib as plt
from scipy.stats import norm

# ファイルをロード
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']
dataframes = [pd.read_csv(f) for f in files]

# 'MPH.csv'を基準として正規化
MPH_df = dataframes[0]
normalized_dfs = [df / MPH_df * 100 for df in dataframes]

# CNS活性のデータを取得 (ここでは仮に'CNS_activity'というカラム名を使用)
# カラム名は各データセットに合わせて変更してください
CNS_activity_data = [df['CNS_activity'] for df in normalized_dfs]

# データを結合
combined_data = np.concatenate(CNS_activity_data)

# データの平均と標準偏差を計算
mu, std = norm.fit(combined_data)

# ヒストグラムと正規分布曲線をプロット
plt.hist(combined_data, bins=30, density=True, alpha=0.6, color='g')

xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', bandwidth=2)
title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)

plt.show()
