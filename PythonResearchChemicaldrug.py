import os
import tkinter as tk
from tkinter import messagebox
from rdkit import Chem
from selenium import webdriver
import pandas as pd
import time

def iupac_to_smiles(iupac):
    mol = Chem.MolFromSmiles(iupac)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def get_adme_data(smiles):
    # SwissADMEからADMEデータを取得する
    # 1.5分待機する
    time.sleep(90)
    # ランダムにプロキシを選択する
    proxies = [ 
         { 'http' : '43.157.66.170:8080' }
        ,{ 'http' : '139.162.78.109:8080' }
        ,{ 'https' : '43.153.222.203:8080' }
        ,{ 'https' : '43.138.216.160:8080' }
        ,{ 'https' : '43.135.182.214:8080' }	

    ]
    proxy = random.choice(proxies)
    # プロキシを使用してSwissADMEにアクセスする
    r = requests.get('http://www.swissadme.ch/', params={'proxy': proxy, 'SMILES': smiles})
    # 結果を取得するまで30秒待機
    time.sleep(30)
adme_data ='SwissADME.csv'
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
    ec50_data = get_ec50_data(adme_data)
    result_text.configure(state='normal')
    result_text.delete(1.0, tk.END)
    result_text.insert(tk.END, {SMILES: {smiles} \Ki:\DAT: {ki_data["DAT"]}\NAT: {ki_data["NAT"]}\SERT: {ki_data["SERT"]}\EC50:\DAT: {ec50_data["DAT"]}\NAT: {ec50_data["NAT"]}\S\NAT: {ki_data["NAT"]}\SERT: {ki_data["SERT"]}EC50:DAT: {ec50_data["DAT"]NAT: {ec50_data["NAT"]}S_data["NAT"]}SERT: {ki_data["SERT"]}EC50:DAT: {ec50_data["DAT"]}NAT: {ec50_data["NAT"]}S
