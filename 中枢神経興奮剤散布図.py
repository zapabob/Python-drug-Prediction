import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、MWとXLOGPのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            data.append((df['MW'], df['XLOGP']))
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, (x, y) in enumerate(data):
        plt.scatter(x, y, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        x = x.values.reshape(-1, 1)
        reg = LinearRegression().fit(x, y)
        plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('MW')
    plt.ylabel('XLOGP')
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
