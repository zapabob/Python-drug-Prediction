import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

# CSVファイルからデータを読み込み
file_names = ['MPH.csv','PR.csv','4MMC.csv','4mar.csv','A.csv','Cocaine.csv','MA.csv']
dataframes = [pd.read_csv(file) for file in file_names]

# データフレームを結合
df = pd.concat(dataframes)

# SMILES（化学構造を表す文字列）から分子指紋を生成
df['mol'] = df['smiles'].apply(Chem.MolFromSmiles)
df['fp'] = df['mol'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024))

# 分子指紋を特徴ベクトルに変換
X = np.array(list(df['fp']))
y = df['property'].values  # 'property' は予測したい物性を表すカラムです

# データを訓練セットとテストセットに分割
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 線形回帰モデルの訓練
model = LinearRegression()
model.fit(X_train, y_train)

# モデルの評価
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse}")
