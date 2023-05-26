import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats

# ファイルのリスト
files = ['MPH.csv','PP.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、MPHに対してMW活性を正規化
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            if file == '4mar.csv':
                mph_MW = df['MW']  # 4marのMWを取得
            else:
                df['Normalized MW'] = df['MW'] / mph_MW
                data.append(df['Normalized MW'])
    return data

# 全てのデータを一つのリストに結合
def combine_data(data):
    return np.concatenate(data)

# ヒストグラムを描き、正規分布をフィット
def plot_distribution(all_data):
    plt.hist(all_data, bins=30, density=True, alpha=0.6, color='g')

    mu, std = scipy.stats.norm.fit(all_data)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = scipy.stats.norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)

    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)
    plt.show()

def main():
    data = process_files(files)
    all_data = combine_data(data)
    plot_distribution(all_data)

if __name__ == "__main__":
    main()
