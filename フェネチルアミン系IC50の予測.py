import pandas as pd
import numpy as np
import time
import requests
import tkinter as tk
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# Prepare and preprocess data
data = {
    "Compound": ["MDMA", "Methylone", "BDB", "MBDB", "MDPV", "2C-I", "2C-E", "2C-C", "TMA", "TMA-2", "TMA-6", "PMMA", "4FMP", "MAP"],
    "DA IC50": [1.4e-6, 2.9e-6, 7.9e-6, 6.3e-6, 4.3e-8, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.4e-5, 7.7e-7, 3.7e-7],
}
df = pd.DataFrame(data)
smiles_dict = {}
for compound in df["Compound"]:
    time.sleep(2)
    try:
        res = requests.get(f"https://api.psychonautwiki.org/{compound}")
        res_json = res.json()
        smiles = res_json["data"]["relationships"]["versions"]["data"][0]["attributes"]["smiles"]
        smiles_dict[compound] = smiles
    except:
        pass
df["SMILES"] = df["Compound"].map(smiles_dict)
df = df.dropna(subset=["SMILES"])
molecules = [Chem.MolFromSmiles(smiles) for smiles in df["SMILES"]]
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in molecules]
fp_df = pd.DataFrame([list(fp) for fp in fingerprints])

# Train model
fp_train, fp_test, ic50_train, ic50_test = train_test_split(fp_df, df["DA IC50"], test_size=0.2, random_state=42)
model = RandomForestRegressor()
model.fit(fp_train, ic50_train)

# Create GUI
def predict_ic50():
    iupac = iupac_entry.get()
    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(iupac))
    molecule = Chem.MolFromSmiles(smiles)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, 2, nBits=1024)
    fp = pd.DataFrame([list(fingerprint)])
    ic50 = model.predict(fp)
    ic50_text.set(f"Predicted IC50: {ic50[0]}")

root = tk.Tk()
root.title("IC50 Predictor") 
root.geometry("300x200")

ic50_text = tk.StringVar()

iupac_label = tk.Label(root, text="Enter IUPAC:") 
iupac_label.pack()
iupac_entry = tk.Entry(root)
iupac_entry.pack()

predict_button = tk.Button(root, text="Predict IC50", command=predict_ic50) 
predict_button.pack()

ic50_label = tk.Label(root, textvariable=ic50_text) 
ic50_label.pack()

root.mainloop()