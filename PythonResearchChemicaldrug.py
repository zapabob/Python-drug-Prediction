import os
import tkinter as tk
from tkinter import messagebox
from rdkit import Chem
from selenium import webdriver
import pandas as pd

def iupac_to_smiles(iupac):
    mol = Chem.MolFromSmiles(iupac)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def get_adme_data(smiles):
    # SwissADMEからADMEデータを取得する
    driver = webdriver.Chrome()
    driver.get('http://www.swissadme.ch/')
    driver.find_element_by_name('SMILES').send_keys(smiles)
    driver.find_element_by_name('qsar').click()
    driver.find_element_by_link_text('Export QSAR results as CSV').click()
    driver.quit()

    # CSVファイルからADMEデータを読み込む
    adme_data = pd.read_csv('SwissADME.csv', delimiter=';')
    return adme_data

def get_ki_data(adme_data):
    # 各モノアミン受容体への半数阻害効果定数、半数効果濃度定数を取得する
    ki_data = {}
    ki_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'Ki (nM)'].iloc[0]
    ki_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'Ki (nM)'].iloc[0]
    ki_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'Ki (nM)'].iloc[0]
    return ki_data

def get_ec50_data(adme_data):
    # 各モノアミン受容体への半数効果濃度定数を取得する
    ec50_data = {}
    ec50_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'EC50 (nM)'].iloc[0]
    ec50_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'EC50 (nM)'].iloc[0]
    ec50_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'EC50 (nM)'].iloc[0]
    return ec50_data

def button_click():
    iupac = iupac_entry.get()
    smiles = iupac_to_smiles(iupac)
    adme_data = get_adme_data(smiles)
    ki_data = get_ki_data(adme_data)
    result_text.configure(state='normal')
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, f'SMILES: {smiles}\n\nDAT: {ki_data["DAT"]}\nNAT: {ki_data["NAT"]}\nSERT: {ki_data["SERT"]}')
    result_text.configure(state='disabled')
    adme_data.to_csv(os.path.expanduser("~/Desktop/My.csv"), index=False)

# GUIアプリケーションのウィンドウを作成する
root = tk.Tk()
root.title('ADME Data')

# IUPAC名を入力するためのテキストボックスを作成し、ウィンドウ上に配置する
iupac_entry = tk.Entry(root, width=40)
iupac_entry.pack(pady=10)

# ボタンを作成し、ウィンドウ上に配置する
button = tk.Button(root, text='Get ADME Data', command=button_click)
button.pack(pady=10)

# 結果を表示するためのテキストボックスを作成し、ウィンドウ上に配置する
result_text = tk.Text(root, width=50, height=10, state='disabled')
result_text.pack(pady=10)

# GUIアプリケーションを実行する
root.mainloop()
