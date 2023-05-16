import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# ファイルのリスト
files = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']

# 各ファイルからデータを読み込み、分子量あたりのMLogPを計算
def process_files(files):
    MW = []
    MLogP_values = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            MW.append(df['MW'])
            MLogP_values.append(df['MLogP'])
    return MW, MLogP_values

# データを散布図で表示し、回帰直線を描く
def plot_data(MW, MLogP_values):
    plt.figure(figsize=(10, 6))
    for i, (w, MLogP) in enumerate(zip(MW, MLogP_values)):
        plt.scatter(w, MLogP, label=files[i].split('.')[0])  # 各CSVファイルのデータを散布図で表示
        reg = LinearRegression().fit(np.array(w).reshape(-1, 1), np.array(MLogP).reshape(-1, 1))
        x = np.linspace(min(w), max(w), 100).reshape(-1, 1)
        plt.plot(x, reg.predict(x), label=f'Regression line {files[i].split(".")[0]}')
    plt.xlabel('MolWeight')
    plt.ylabel('MLogP')
    plt.legend()
    plt.show()

def main():
    MW, MLogP_values = process_files(files)
    plot_data(MW, MLogP_values)

if __name__ == "__main__":
    main()
