from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense
import tkinter as tk
from rdkit.Chem import MolFromSmiles, MolToSmiles
import pandas as pd
from sklearn.metrics import r2_score
# ChEMBLデータベースからデータを取得
target = new_client.target
target_query = target.search('CHEMBL238')
targets = pd.DataFrame.from_dict(target_query)

selected_target = targets.target_chembl_id[0]
activities = new_client.activity
res = activities.filter(target_chembl_id=selected_target).filter(standard_type='IC50')

df = pd.DataFrame.from_dict(res)

# 分子記述子を計算
def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    fp_arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, fp_arr)

    descriptors = np.concatenate(([mol_weight, logp, num_h_donors, num_h_acceptors], fp_arr))

    return descriptors

# データの準備
data = []
labels = []
for _, act in df.iterrows():
    if act['standard_value'] is not None and act['canonical_smiles'] is not None:
        descriptors = compute_descriptors(act['canonical_smiles'])
        if descriptors is not None:
            data.append(descriptors)
            labels.append(float(act['standard_value']))

data = np.array(data)
labels = np.array(labels)

# データの分割
train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.2, random_state=42)

# モデルの構築
model = Sequential([
    Dense(1024, activation='relu', input_shape=(train_data.shape[1],)),
    Dense(512, activation='relu'),
    Dense(1)
])

# モデルのコンパイル
model.compile(optimizer='adam', loss='mean_squared_error')

# モデルの訓練
model.fit(train_data, train_labels, epochs=1000, batch_size=32, validation_split=0.2)

# モデルの評価
test_loss = model.evaluate(test_data, test_labels)
print(f"Test loss: {test_loss}")




# モデルの予測を行う関数
def predict_ic50(smiles):
    # SMILESを分子記述子に変換
    descriptors = compute_descriptors(smiles)
    # モデルによる予測
    predicted_ic50 = model.predict(np.array([descriptors]))
    # IC50の予測値を-pIC50として返す
    return -np.log10(predicted_ic50)

# ボタンが押されたときの処理
def on_button_press():
    smiles = iupac_entry.get()
    try:
        predicted_ic50 = predict_ic50(smiles)
        result_label.config(text=f"Predicted -pIC50: {predicted_ic50}")
    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")
# モデルの予測を取得
predictions = model.predict(test_data)

# R^2スコアを計算
r2 = r2_score(test_labels, predictions)

print(f"R^2 score: {r2}")

    # GUIの作成
root = tk.Tk()
root.title("IC50 Predictor")

    # 入力フィールド
iupac_entry = tk.Entry(root)
iupac_entry.pack()

    # 結果表示ラベル
result_label = tk.Label(root)
result_label.pack()
    # ボタン
predict_button = tk.Button(root, text="Predict -pIC50", command=on_button_press)
predict_button.pack()

root.mainloop()