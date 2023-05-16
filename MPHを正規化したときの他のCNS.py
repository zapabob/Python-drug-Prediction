import numpy as np
import pandas as pd
import matplotlib.pyplot as plt  # <- "matplotlib.pyplot"をインポート
import scipy.stats

# CSVファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、MPHに対してCNS活性を正規化
data = []
for file in files:
    df = pd.read_csv(file)
    if file == 'MPH.csv':
        mph_cns = df['CNS']  # MPHのCNS活性値を取得
    else:
        df['Normalized CNS'] = df['CNS'] / mph_cns  # <- CNSを文字列として指定
        data.append(df['Normalized CNS'])

# 全てのデータを一つのリストに結合
all_data = np.concatenate(data)

# ヒストグラムを描く
plt.hist(all_data, bins=30, density=True, alpha=0.6, color='g')  # <- "plt.hist"を使用

# 正規分布のフィット
mu, std = scipy.stats.norm.fit(all_data)
xmin, xmax = plt.xlim()  # <- "plt.xlim"を使用
x = np.linspace(xmin, xmax, 100)
p = scipy.stats.norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)  # <- "linewidth"を使用

# タイトルを設定して描画
title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)
plt.show()
