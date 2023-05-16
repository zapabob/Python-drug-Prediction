import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from rdkit import Chem
from rdkit.Chem import Descriptors

# Step 1: データの準備
data1 = pd.read_csv("MA.csv")
data2 = pd.read_csv("MPH.csv")
data3 = pd.read_csv("Cocaine.csv")

# 複数のcsvファイルを結合
data = pd.concat([data1, data2, data3])

# 欠損値処理
data = data.fillna(data.mean())

# Step 2: 特徴量の抽出
# IUPAC名から特徴量を抽出するためにRDKitを使用
def get_descriptors(iupac_name):
    mol = Chem.MolFromSmiles(iupac_name)
    return Descriptors.MolWt(mol), Descriptors.MolLogP(mol)  # ここでは分子量とLogPを特徴量として使用

data["MolWt"], data["MolLogP"] = zip(*data["IUPAC_name"].map(get_descriptors))

X = data[["MolWt", "MolLogP"]]  # 特徴量
y = data["CNS_activity"]  # 目的変数

# Step 3: モデルの訓練
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = LinearRegression()
model.fit(X_train, y_train)

# Step 4: モデルの評価
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse}")
