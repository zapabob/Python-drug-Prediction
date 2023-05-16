import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、MWあたりのXLOGPのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            XLOGP_per_MW = df['XLOGP'] * df['MW']  # MW*XLOGPを計算
            data.append(XLOGP_per_MW)
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    x = np.array(range(len(data))).reshape(-1, 1)
    for i, d in enumerate(data):
        plt.scatter([i+1]*len(d), d, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(x, d.values.reshape(-1, 1))
        plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('File')
    plt.ylabel('XLOGP * MW')
    plt.xticks(range(1, len(files)+1), [f.split('.')[0] for f in files])
    plt.xlim([0, len(files)+1])  # x軸の表示範囲を設定
    plt.ylim([min(min(data)), max(max(data))])  # y軸の表示範囲を設定
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()


# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, (logMW, XLOGP) in enumerate(data):
        plt.scatter(logMW, XLOGP, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(logMW.values.reshape(-1, 1), XLOGP.values.reshape(-1, 1))
        predicted_values = reg.predict(logMW.values.reshape(-1, 1))
        plt.plot(logMW, predicted_values, color='red')
    plt.xlabel('LogMW')
    plt.ylabel('XLOGP')
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
