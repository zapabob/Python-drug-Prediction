import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、XLOGPとMWのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            data.append(df[['MW', 'XLOGP']])
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, df in enumerate(data):
        x = df['MW'].values.reshape(-1, 1)
        y = df['XLOGP'].values.reshape(-1, 1)
        plt.scatter(np.log10(x), np.log10(y), label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(np.log10(x), np.log10(y))
        x_line = np.linspace(np.min(np.log10(x)), np.max(np.log10(x)), 100).reshape(-1, 1)
        plt.plot(x_line, reg.predict(x_line), color='red')
    plt.xlabel('Log10(MW)')
    plt.ylabel('Log10(XLOGP)')
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()


def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
