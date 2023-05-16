import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、XLOGPとLOG(MW)のデータを取得
def process_files(files):
    data_x = []
    data_y = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            data_x.append(np.log(df['MW']))  # LOG(MW)を計算
            data_y.append(df['XLOGP'])
    return data_x, data_y

# データを散布図で表示し、回帰直線を描く
def plot_data(data_x, data_y):
    plt.figure(figsize=(10, 6))
    for i in range(len(data_x)):
        x = data_x[i].values.reshape(-1, 1)
        y = data_y[i].values.reshape(-1, 1)
        plt.scatter(x, y, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(x, y)
        plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('Log(MW)')
    plt.ylabel('XLOGP')
    plt.xlim([np.min([np.min(dx) for dx in data_x]), np.max([np.max(dx) for dx in data_x])])  # x軸の表示範囲を設定
    plt.ylim([np.min([np.min(dy) for dy in data_y]), np.max([np.max(dy) for dy in data_y])])  # y軸の表示範囲を設定
    plt.legend()
    plt.show()

def main():
    data_x, data_y = process_files(files)
    plot_data(data_x, data_y)

if __name__ == "__main__":
    main()
