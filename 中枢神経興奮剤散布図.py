import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv']

# 各ファイルからデータを読み込み、LOGSとXLOGP3のデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            data.append(df[['LOGS', 'XLOGP3']])
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, df in enumerate(data):
        x = df['LOGS'].values.reshape(-1, 1)
        y = df['XLOGP3'].values.reshape(-1, 1)
        plt.scatter(x, y, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        if len(x) > 1:  # 回帰直線を描くためには少なくとも2つのデータ点が必要です
            reg = LinearRegression().fit(x, y)
            plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('LOGS')
    plt.ylabel('XLOGP3')
    plt.xlim([np.min([df['LOGS'].min() for df in data]), np.max([df['LOGS'].max() for df in data])])  # x軸の表示範囲を設定
    plt.ylim([np.min([df['XLOGP3'].min() for df in data]), np.max([df['XLOGP3'].max() for df in data])])  # y軸の表示範囲を設定
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
