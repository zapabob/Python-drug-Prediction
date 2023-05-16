import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、logMWとXLOGPのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            logMW = np.log(df['MW'])  # logMWを計算
            XLOGP = df['XLOGP']  # XLOGPを取得
            data.append((logMW, XLOGP))
    return data

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
