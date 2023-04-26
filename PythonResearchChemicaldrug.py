#!/bin/env python3
import random
import pandas as pd
from classes import MainWindow

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


if __name__ == "__main__":
    main()