#!/bin/env python3

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions
import re
import requests
import io
import sys

class SwissADME:
    def __init__(self, smiles, proxy=None):
        # Chromeのヘッドレスモードを有効にするオプションを設定します
        chrome_options = webdriver.ChromeOptions()
        chrome_options.add_argument("--headless")

        # プロキシ設定
        self.proxy = None
        if proxy:
            m = re.match(r"https*://[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+:[0-9]+", proxy)
            if m:
                self.proxy = proxy
                print(f"proxy setting is enabled. ({self.proxy})")
                chrome_options.add_argument(f"--proxy-server={self.proxy}")
            else:
                raise Exception("Proxy server is not found. Please confirm your setting.")

        # ChromeDriverを初期化し、ブラウザを起動します
        driver = webdriver.Chrome(options=chrome_options)

        # SwissADMEのWebサイトにアクセスします
        driver.get("http://www.swissadme.ch/index.php")

        # テキストボックスに文字列を入力します
        text_box = driver.find_element(By.ID, "smiles")
        text_box.send_keys(smiles)

        # ボタンを押下します
        submit_button = driver.find_element(By.ID, "submitButton")
        submit_button.click()

        # ページが更新されるのを待ち、CSVダウンロードリンクが表示されるまで待ちます
        WebDriverWait(driver, 30).until(
            expected_conditions.presence_of_element_located((By.CLASS_NAME, "expendable"))
        )

        # HTMLソースを取得
        html = driver.page_source.replace("\n", "")

        # CSVファイルのダウンロードリンクを取得
        url_csv = ""
        m = re.match(r'.*(results/.*/swissadme\.csv).*', html)
        if m:
            url_csv = "http://www.swissadme.ch/" + m.group(1)
        else:
            raise Exception("There is no CSV link.")

        # CSVファイルをダウンロード
        response = requests.get(url_csv, proxies=self.proxy)
        if response.status_code != 200:
            raise Exception("Download Error!")

        # ダウンロードしたCSVファイルをDataFrameとして読み込みます
        csv_strio = io.StringIO(response.text)
        self.data_frame = pd.read_csv(csv_strio)

        # ウェブドライバを閉じてリソースを解放します
        driver.quit()

    def get(self) -> pd.DataFrame:
        return self.data_frame

def main():
    smiles: str = ""

    # コマンドライン引数チェック
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <SMILES description>")
        print('デフォルト動作として、Methamphetamine: ”CNC(C)Cc1ccccc1”の情報を取得します')
        smiles = "CNC(C)Cc1ccccc1"
    else:
        smiles = sys.argv[1]

    # ADMEデータを取得
    adme = SwissADME(smiles)
    df: pd.DataFrame = adme.get()

    # 表示
    for _, row in df.iterrows():
        print(f"{row.name}: {row}")

if __name__ == "__main__":
    main()
