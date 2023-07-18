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

# ChEMBLデータベースからデータを取得
target = new_client.target
target_query = target.search('CHEMBL238')
targets = pd.DataFrame.from_dict(target_query)
targets

selected_target = targets.target_chembl_id
selected_target
activities = new_client.activity
res = activities.filter(target_chembl_id='CHEMBL238').filter(standard_type='IC50')
res
df = pd.DataFrame.from_dict(res)
df
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
for act in res:
    if act['standard_value'] is not None and act['canonical_smiles'] is not None:
        descriptors = compute_descriptors(act['canonical_smiles'])
        if descriptors is not None:
            data.append(descriptors)
            labels.append(float(act['standard_value']))

data = np.array(data)
labels = np.array(labels)

# データの分割
train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.2, random_state=64)

# モデルの構築
from keras.layers import Dropout

# モデルの構築
model = Sequential([
    Dense(1024, activation='relu', input_shape=(train_data.shape[1],)),
    Dropout(0.5),
    Dense(512, activation='relu'),
    Dropout(0.5),
    Dense(1)
])


# モデルのコンパイル
model.compile(optimizer='adam', loss='mean_squared_error')

# モデルの訓練
model.fit(train_data, train_labels, epochs=25, batch_size=64, validation_split=0.2)

# モデルの評価
test_loss = model.evaluate(test_data, test_labels)
print(f"Test loss: {test_loss}")
# モデルの予測
test_predictions = model.predict(test_data)


# モデルの予測を行う関数
def predict_ic50(iupac_name):
    # IUPAC名をSMILESに変換
    smiles = MolToSmiles(MolFromSmiles(iupac_name))
    # SMILESを分子記述子に変換
    descriptors = compute_descriptors(smiles)
    # モデルによる予測
    predicted_ic50 = model.predict(np.array([descriptors]))
    # IC50が1000を超える場合はN/Aを返す
    if predicted_ic50 > np.log10(10^-3):
        return "N/A" # IC50 >1000
    else:
        return -np.log10(predicted_ic50)

from rdkit import Chem
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt

compounds = {
    "Cocaine": 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)C3=CC=CC=C3)C(=O)OC',
    "Methamphetamine": 'CNC(C)Cc1ccccc1',
    "MDPV": 'O=C(C(CCC)N1CCCC1)C2=CC=C3C(OCO3)=C2',
    "GBR 12909":'Fc1ccc(cc1)C(OCCN2CCN(CC2)CCCc3ccccc3)c4ccc(F)cc4',
    "a-PHP":'C1(=CC=CC=C1)C(C(CCCC)N2CCCC2)=O',
    "a-PVP":'CCCC(C(C1=CC=CC=C1)=O)N2CCCC2',
    "4-mar":'CC1C(C2=CC=CC=C2)OC(N)=N1',
    "4Br4MAR":'Brc1ccc(cc1)C1OC(N)=NC1C',
    "4.4DMAR":'CC(N=C(N)O1)C1C2=CC=C(C)C=C2',
    "MPH":'COC(=O)C(c1ccccc1)C1CCCCN1',
    "MDMA":'CC(NC)CC1=CC=C(OCO2)C2=C1',
    "AMP":'C[C@@H](Cc1ccccc1)N',
}

logP_values = []
predicted_pic50_values = []

for name, smiles in compounds.items():
    mol = Chem.MolFromSmiles(smiles)
    logP = Descriptors.MolLogP(mol)
    predicted_pic50 = abs(predict_ic50(smiles))
    
    logP_values.append(logP)
    predicted_pic50_values.append(predicted_pic50)
    
    print(f"LogP and predicted pIC50 for {name}: {logP}, {predicted_pic50}")

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(logP_values, predicted_pic50_values)
for i, txt in enumerate(compounds.keys()):
    plt.annotate(txt, (logP_values[i], predicted_pic50_values[i]))
plt.xlabel('LogP')
plt.ylabel('Predicted pIC50')
plt.title('Scatter plot of predicted pIC50 against LogP')
plt.show()
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np

# Convert the lists to numpy arrays and reshape them
logP_values_np = np.array(logP_values).reshape(-1, 1)
predicted_pic50_values_np = np.array(predicted_pic50_values).reshape(-1, 1)

# Fit the linear regression model
model = LinearRegression()
model.fit(logP_values_np, predicted_pic50_values_np)

# Predict the pIC50 values
predicted_pic50_values_pred = model.predict(logP_values_np)

# Compute the R^2 score
r2 = r2_score(predicted_pic50_values_np, predicted_pic50_values_pred)

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(logP_values, predicted_pic50_values)
plt.plot(logP_values, predicted_pic50_values_pred, color='red')  # Add the regression line
for i, txt in enumerate(compounds.keys()):
    plt.annotate(txt, (logP_values[i], predicted_pic50_values[i]))
plt.xlabel('LogP')
plt.ylabel('Predicted pIC50')
plt.title(f'Scatter plot of predicted pIC50 against LogP with regression line (R^2 = {r2:.2f})')
plt.show()

# Print the equation of the line
print(f"Equation of the line: pIC50 = {model.intercept_[0]:.2f} + {model.coef_[0][0]:.2f} * LogP")


# GUIの作成
root = tk.Tk()
root.title("IC50 Predictor")


entry = tk.Entry(root)
entry.pack()


def reset():
    # clear any stateful parts of the application
    entry.delete(0, 'end')
    result_label.config(text="")
    # clear the input field
    iupac_entry.delete(0, 'end')
    # you might also need to reset your model or any other stateful parts of your application here


reset_button = tk.Button(root, text="Reset", command=reset)
reset_button.pack()



def copy():
    root.clipboard_clear()
    root.clipboard_append(entry.get())

def paste():
    entry.insert(0, root.clipboard_get())

copy_button = tk.Button(root, text="Copy", command=copy)
copy_button.pack()

paste_button = tk.Button(root, text="Paste", command=paste)
paste_button.pack()


# 入力フィールド
iupac_entry = tk.Entry(root)
iupac_entry.pack()

# 結果表示ラベル
result_label = tk.Label(root)
result_label.pack()

# ボタンが押されたときの処理
def on_button_press():
    iupac_name = iupac_entry.get()
    try:
        predicted_ic50 = -predict_ic50(iupac_name)
        result_label.config(text=f"Predicted IC50: {predicted_ic50}")
    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")

# ボタン
predict_button = tk.Button(root, text="Predict pIC50", command=on_button_press)
predict_button.pack()

root.mainloop()