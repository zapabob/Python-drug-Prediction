from rdkit.Chem import Descriptors
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from rdkit import Chem

# List of files
files = ['MPH.csv','MA.csv','A.csv']

def process_files(files):
    data = []
    for file in files:
        if os.path.isfile(file):
            df = pd.read_csv(file)
            df['Molecule'] = df['SMILES'].apply(Chem.MolFromSmiles)
            df['MW'] = df['Molecule'].apply(Descriptors.MolWt)
            df['TPSA'] = df['Molecule'].apply(Descriptors.TPSA)
            data.append(df[['TPSA','MW']])
    return data

def plot_data(data):
    plt.figure(figsize=(10, 6))
    for i, df in enumerate(data):
        x = df['TPSA'].values.reshape(-1, 1)
        y = df['MW'].values.reshape(-1, 1)
        plt.scatter(x, y, label=files[i].split('.')[0])  # Display data from each CSV file as a scatter plot
        if len(x) > 1:  # At least two data points are required to draw a regression line
            reg = LinearRegression().fit(x, y)
            plt.plot(x, reg.predict(x), color='red')
    plt.xlabel('TPSA')
    plt.ylabel('MW')
    plt.xlim([np.min([df['TPSA'].min() for df in data]), np.max([df['TPSA'].max() for df in data])])  # Set display range for x-axis
    plt.ylim([np.min([df['MW'].min() for df in data]), np.max([df['MW'].max() for df in data])])  # Set display range for y-axis
    plt.legend()
    plt.show()

def main():
    data = process_files(files)
    plot_data(data)

if __name__ == "__main__":
    main()
