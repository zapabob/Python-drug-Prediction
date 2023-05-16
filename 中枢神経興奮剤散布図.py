import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、MWとXLogPのデータを取得
def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            data.append(df[['MW', 'XLogP']])
    return data

# データを散布図で表示し、回帰直線を描く
def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, d in enumerate(data):
        plt.scatter(d['MW'], d['XLogP'], label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        x = d['MW'].values.reshape(-1, 1)
        y = d['XLogP'].values.reshape(-1, 1)
        reg = LinearRegression().fit(x, y)
        plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('MW')
    plt.ylabel('XLogP')
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
