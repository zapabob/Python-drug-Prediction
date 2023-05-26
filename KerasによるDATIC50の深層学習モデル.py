from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import numpy as np
import tensorflow as tf
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense

# ChEMBLデータベースからデータを取得
target = new_client.target
target_query = target.search('CHEMBL238')
targets = target_query.all()

activities = new_client.activity
res = activities.filter(target_chembl_id=targets[0]['target_chembl_id']).filter(standard_type="IC50")

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
    

    descriptors = np.concatenate(([mol_weight, logp, num_h_donors, num_h_acceptors], fp_arr))

    return descriptors

# データの準備
data = []
labels = []
for act in res:
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
model.fit(train_data, train_labels, epochs=1000, batch_size=3200, validation_split=0.2)

# モデルの評価
test_loss = model.evaluate(test_data, test_labels)
print(f"Test loss: {test_loss}")
import tkinter as tk
from tkinter import messagebox

def predict_ic50():
    smiles = entry.get()
    descriptors = compute_descriptors(smiles)
    if descriptors is not None:
        descriptors = np.array([descriptors])
        prediction = model.predict(descriptors)
        messagebox.showinfo("Prediction", f"The predicted IC50 value is {prediction[0][0]}")
    else:
        messagebox.showerror("Error", "Failed to compute descriptors for the given SMILES string.")

import tkinter as tk
from rdkit import Chem

def predict_ic50():
    # IUPAC名を取得
    iupac_name = iupac_entry.get()

    # IUPAC名をSMILESに変換
    mol = Chem.MolFromSmiles(iupac_name)
    smiles = Chem.MolToSmiles(mol)

    # SMILESから分子記述子を計算
    descriptors = compute_descriptors(smiles)

    # モデルを使用してIC50の値を予測
    ic50 = model.predict(np.array([descriptors]))

    # IC50の値を表示
    result_label.config(text=f"Predicted IC50: {ic50[0][0]}")

# GUIを作成
root = tk.Tk()

# IUPAC名を入力するためのテキストボックスを作成
iupac_entry = tk.Entry(root)
iupac_entry.pack()

# IC50の値を予測するためのボタンを作成
predict_button = tk.Button(root, text="Predict IC50", command=predict_ic50)
predict_button.pack()

# 予測結果を表示するためのラベルを作成
result_label = tk.Label(root, text="")
result_label.pack()

# GUIを実行
root.mainloop()

