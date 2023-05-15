#!/bin/env python3
import pandas as pd
import tkinter as tk
import classes as cl
import threading
import time

def get_adme_data(smiles: str) -> pd.DataFrame:
    # SwissADMEからADMEデータを取得する
    adme = cl.SwissADME(smiles)
    return adme.get()

def is_prime(number):
    if number < 2:
        return False
    for i in range(2, number):
        if number % i == 0:
            return False
    return True

def roll_dice():
    return random.randint(1, 6)

while True:
    dice_result = roll_dice()
    print(f"サイコロの結果: {dice_result}")

    if is_prime(dice_result):
        break

    wait_time = dice_result * 10
    print(f"{wait_time}秒待機します...")
    time.sleep(wait_time)

print("アクセスを再開します。")
def get_ic50_data(adme_data: pd.DataFrame) -> dict[str, str]:
    # 各モノアミン受容体への半数阻害効果濃度を取得する
    IC50_data = {}
    IC50_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'Ki (nM)'].iloc[0]
    IC50_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'Ki (nM)'].iloc[0]
    IC50_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'Ki (nM)'].iloc[0]
    return IC50_data


def get_ec50_data(adme_data: pd.DataFrame) -> dict[str, str]:
    # 各モノアミン受容体への半数効果濃度を取得する
    ec50_data = {}
    ec50_data['DAT'] = adme_data.loc[adme_data['Target'] == 'DAT', 'EC50 (nM)'].iloc[0]
    ec50_data['NAT'] = adme_data.loc[adme_data['Target'] == 'NET', 'EC50 (nM)'].iloc[0]
    ec50_data['SERT'] = adme_data.loc[adme_data['Target'] == 'SERT', 'EC50 (nM)'].iloc[0]
    return ec50_data


# the main function
# This function has no arguments and returns no values.
def main() -> None:
    # make the main window.
    MainWindow("ADME Data")


    # IUPAC->SMILES->ADME Data
    iupac: str = iupac_entry.get()
    smiles: str = cl.SubstanceName(iupac=iupac).get()
    adme_data: pd.DataFrame = get_adme_data(smiles)

    # ki_data = get_ic50_data(adme_data)
    # ec50_data = get_ec50_data(adme_data)

    # 結果の表示
    result_text.configure(state='normal')
    result_text.delete(1.0, tk.END)
    # result_text.insert(tk.END, f'SMILES: {smiles}\n\nDAT:\n  Ki: {ki_data["DAT"]} nM\n  EC50: {ec50_data["DAT"]} nM\n\nNAT:\n  Ki: {ki_data["NAT"]} nM\n  EC50: {ec50_data["NAT"]} nM\n\nSERT:\n  Ki: {ki_data["SERT"]} nM\n  EC50: {ec50_data["SERT"]} nM')
    text: str = ""
    for _, item in adme_data.iterrows():
        text = str(item)
    result_text.insert(tk.END, text)
    result_text.configure(state='disabled')

    # CSV出力
    # adme_data.to_csv("adme.csv", encoding="utf_8_sig")

    # ボタンを有効化する
    thread = threading.Thread(target=activate_button, args=[10])
    thread.start()

def activate_button(dulation: int) -> None:
    # dulationで指定された秒数でボタンを有効化する
    for i in range(dulation, 0, -1):
        button.configure(text=f"{i}")
        time.sleep(1)

    button.configure(text="Get ADME Data")
    button.configure(state="normal")

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
