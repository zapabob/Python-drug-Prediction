from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MolFromSmiles, MolToSmiles
import numpy as np
from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense, Dropout
import tkinter as tk
import pandas as pd
import matplotlib.pyplot as plt

# ChEMBLデータベースからデータを取得
target = new_client.target

# 複数のターゲットを個別に検索
target_ids = ['CHEMBL218']
all_res = []

for target_id in target_ids:
    target_query = target.search(target_id)
    selected_target = target_query[0]['target_chembl_id']
    activities = new_client.activity
    res = activities.filter(target_chembl_id=selected_target).filter(standard_type='Ki')
    all_res.extend(res)

df = pd.DataFrame.from_dict(all_res)

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    fp_arr = np.zeros((1,))
    AllChem.DataStructs.ConvertToNumpyArray(fp, fp_arr)

    descriptors = np.concatenate(([mol_weight, logp, num_h_donors, num_h_acceptors], fp_arr))
    return descriptors
data = []
labels = []
for act in all_res:
    if act['standard_value'] is not None and act['canonical_smiles'] is not None:
        descriptors = compute_descriptors(act['canonical_smiles'])
        if descriptors is not None:
            data.append(descriptors)
            labels.append(float(act['standard_value']))

data = np.array(data)
labels = np.array(labels)
train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.2, random_state=64)
compounds = {
    "Nabilon":'CCCCCCC(C)(C)c1cc(O)c2c(c1)OC(C)(C)C1CCC(=O)CC21',
    "CP55940":'CCCCCCC(C)(C)c1ccc([C@@H]2C[C@H](O)CC[C@H]2CCCO)c(O)c1',
    "D9THC":'CCCCCc1cc(O)c2[C@@H]3C=C(C)CC[C@H]3C(C)(C)Oc2c1',
    "Anandamide":'O=C(NCCO)CCC\C=C/C\C=C/C\C=C/C\C=C/CCCCC',
    "D9THCH":'CCCCCCC1=CC(=C2[C@@H]3C=C(CC[C@H]3C(OC2=C1)(C)C)C)O',
    "JWH018":'CCCCCN1C=C(C(C2=CC=CC3=CC=CC=C32)=O)C4=CC=CC=C41',
    "5-CL-ADB-A":'O=C(N[C@@H](C(C)(C)C)C(OC)=O)C1=NN(CCCC=C)C2=C1C=CC=C2',
}
# モデルの構造を修正
model = Sequential([
    Dense(1024, activation='relu', input_shape=(train_data.shape[1],)),
    Dropout(0.5),
    Dense(512, activation='relu'),
    Dropout(0.5),
    Dense(1)  # 出力層のノード数を1に修正
])

model.compile(optimizer='adam', loss='mean_squared_error')
model.fit(train_data, train_labels, epochs=25, batch_size=64, validation_split=0.2)

test_loss = model.evaluate(test_data, test_labels)
print(f"Test loss: {test_loss}")

def predicted_Ki(smiles):
    descriptors = compute_descriptors(smiles)
    predicted_value = model.predict(np.array([descriptors]))[0][0]
    return np.log10(predicted_value)

logP_values = []
predicted_Ki_values = []

for name, smiles in compounds.items():
    mol = Chem.MolFromSmiles(smiles)
    logP = Descriptors.MolLogP(mol)
    predicted_value = predicted_Ki(smiles)
    
    if predicted_value != "N/A":
        logP_values.append(logP)
        predicted_Ki_values.append(abs(float(predicted_value)))
# compounds辞書の各化合物のKiの予測値を計算
predicted_Ki_values_dict = {name: predicted_Ki(smiles) for name, smiles in compounds.items()}

# 予測されたKi値を出力
for compound, predicted_value in predicted_Ki_values_dict.items():
    print(f"{compound}: Predicted Ki = {predicted_value}")

plt.figure(figsize=(10, 6))
plt.scatter(logP_values, predicted_Ki_values)
for i, (name, smiles) in enumerate(compounds.items()):
    predicted_value = predicted_Ki(smiles)
    if predicted_value != "N/A":
        plt.annotate(name, (logP_values[i], predicted_Ki_values[i]))

plt.xlabel('LogP')
plt.ylabel('Predicted Ki')
plt.title('Scatter plot of predicted Ki against LogP')
plt.yscale('log')  # 縦軸を対数スケールに設定
plt.tight_layout()
plt.show()
root = tk.Tk()
root.title("Ki Predictor")

entry = tk.Entry(root)
entry.pack()

def reset():
    entry.delete(0, 'end')
    result_label.config(text="")
    iupac_entry.delete(0, 'end')

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

iupac_entry = tk.Entry(root)
iupac_entry.pack()

result_label = tk.Label(root)
result_label.pack()

def on_button_press():
    iupac_name = iupac_entry.get()
    try:
        predicted_value = predicted_Ki(iupac_name)
        result_label.config(text=f"Predicted Ki: {predicted_value}")
    except Exception as e:
        result_label.config(text=f"Error: {str(e)}")

predict_button = tk.Button(root, text="Predict ki", command=on_button_press)
predict_button.pack()

root.mainloop()