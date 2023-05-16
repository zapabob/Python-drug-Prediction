import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、分子量あたりのXLoGPを計算
def process_files(files):
    MW = []
    XLoGP_values = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            MW.append(df['MW'])
            XLoGP_values.append(df['XLoGP'])
    return MW, XLoGP_values

# データを散布図で表示し、回帰直線を描く
def plot_data(MW, XLoGP_values):
    plt.figure(figsize=(10, 6))
    for i, (MW, XLoGP) in enumerate(zip(MW, XLoGP_values)):
        plt.scatter(MW, XLoGP, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(np.array(MW).reshape(-1, 1), np.array(XLoGP).reshape(-1, 1))
        x = np.linspace(min(MW), max(MW), 100).reshape(-1, 1)
        plt.plot(x, reg.predict(x), label=f'Regression line {files[i].split(".")[0]}')
    plt.xlabel('MW')
    plt.ylabel('XLoGP')
    plt.legend()
    plt.shoMW()

def main():
    MW, XLoGP_values = process_files(files)
    plot_data(MW, XLoGP_values)

if __name__ == "__main__":
    main()
