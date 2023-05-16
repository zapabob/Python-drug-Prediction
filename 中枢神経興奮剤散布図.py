import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、e(-XLOGP3)とMWのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            df['e(-XLOGP3)'] = np.exp(-df['XLOGP3'])  # e(-XLOGP3)を計算
            data.append(df[['e(-XLOGP3)', 'MW']])
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, df in enumerate(data):
        x = df['MW'].values.reshape(-1, 1)
        y = df['e(-XLOGP3)'].values.reshape(-1, 1)
        plt.scatter(x, y, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(x, y)
        plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('MW')
    plt.ylabel('e(-XLOGP3)')
    plt.xlim([np.min([df['MW'].min() for df in data]), np.max([df['MW'].max() for df in data])])  # x軸の表示範囲を設定
    plt.ylim([np.min([df['e(-XLOGP3)'].min() for df in data]), np.max([df['e(-XLOGP3)'].max() for df in data])])  # y軸の表示範囲を設定
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
