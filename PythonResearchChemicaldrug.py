#!/bin/env python3

import pandas as pd
import tkinter as tk
import cirpy
import asyncio
from classes.SwissADME import SwissADME

def iupac_to_smiles(iupac: str) -> str:
    return cirpy.resolve(iupac, "smiles")

def get_adme_data(smiles: str) -> pd.DataFrame:
    # SwissADMEからADMEデータを取得する
    adme = SwissADME(smiles)
    return adme.get()

def get_ki_data(adme_data: pd.DataFrame) -> dict[str, str]:
    # 各モノアミン受容体への半数阻害効果定数、半数効果濃度定数を取得する
    ki_data = {}
    ki_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'Ki (nM)'].iloc[0]
    ki_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'Ki (nM)'].iloc[0]
    ki_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'Ki (nM)'].iloc[0]
    return ki_data

def get_ec50_data(adme_data: pd.DataFrame) -> dict[str, str]:
    # 各モノアミン受容体への半数効果濃度定数を取得する
    ec50_data = {}
    ec50_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'EC50 (nM)'].iloc[0]
    ec50_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'EC50 (nM)'].iloc[0]
    ec50_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'EC50 (nM)'].iloc[0]
    return ec50_data

def button_click() -> None:
    iupac: str = iupac_entry.get()
    smiles: str = iupac_to_smiles(iupac)
    adme_data: pd.DataFrame = get_adme_data(smiles)
    # ki_data = get_ki_data(adme_data)
    # ec50_data = get_ec50_data(adme_data)
    result_text.configure(state='normal')
    result_text.delete(1.0, tk.END)
    # result_text.insert(tk.END, f'SMILES: {smiles}\n\nDAT:\n  Ki: {ki_data["DAT"]} nM\n  EC50: {ec50_data["DAT"]} nM\n\nNAT:\n  Ki: {ki_data["NAT"]} nM\n  EC50: {ec50_data["NAT"]} nM\n\nSERT:\n  Ki: {ki_data["SERT"]} nM\n  EC50: {ec50_data["SERT"]} nM')
    text: str = ""
    for _, item in adme_data.iterrows():
        text = str(item)
    result_text.insert(tk.END, text)
    result_text.configure(state='disabled')
    adme_data.to_csv("adme.csv", encoding="utf_8_sig")

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
result_text = tk.Text(root, width=80, height=20, state='disabled')
result_text.pack(pady=10)

# GUIアプリケーションを実行する
root.mainloop()

