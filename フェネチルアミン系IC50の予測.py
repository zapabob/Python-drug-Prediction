import pandas as pd
import numpy as np
import time
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from tkinter import Tk, Label, Entry, Button, StringVar

# DataFrame
data = {
    "Compound": ["MDMA", "Methylone", "BDB", "MBDB", "MDPV", "2C-I", "2C-E", "2C-C", "TMA", "TMA-2", "TMA-6", "PMMA", "4FMP", "MAP"],
    "DA IC50": [1.4e-6, 2.9e-6, 7.9e-6, 6.3e-6, 4.3e-8, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.4e-5, 7.7e-7, 3.7e-7],
    # Add other measurements here as necessary
}

df = pd.DataFrame(data)

# Get SMILES strings from the API
smiles_dict = {}
for compound in df["Compound"]:
    time.sleep(2)  # wait for 2 seconds between each API call
    try:
        res = requests.get(f"https://api.psychonautwiki.org/{compound}")
        res_json = res.json()
        smiles = res_json["data"]["relationships"]["versions"]["data"][0]["attributes"]["smiles"]
        smiles_dict[compound] = smiles
    except:
        pass

# Add SMILES strings to the DataFrame
df["SMILES"] = df["Compound"].map(smiles_dict)

# Drop rows with missing SMILES strings
df = df.dropna(subset=["SMILES"])

# Compute fingerprints
molecules = [Chem.MolFromSmiles(smiles) for smiles in df["SMILES"]]
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in molecules]

# Convert fingerprints to a DataFrame
fp_df = pd.DataFrame([list(fp) for fp in fingerprints])

# Split the dataset into a training set and a test set
fp_train, fp_test, ic50_train, ic50_test = train_test_split(fp_df, df["DA IC50"], test_size=0.2, random_state=42)

# Train a model
model = RandomForestRegressor()
model.fit(fp_train, ic50_train)

# Create a GUI for predictions
root = Tk()

input_label = Label(root, text="Enter IUPAC name:")
input_label.pack()

input_text = StringVar()
input_entry = Entry(root, textvariable=input_text)
input_entry.pack()

output_label = Label(root, text="Predicted DA IC50:")
output_label.pack()

output_text = StringVar()
output_entry = Entry(root, textvariable=output_text)
output_entry.pack()

def make_prediction():
    compound = input_text.get()

res = requests.get(f"https://api.psychonautwiki.org/{compound}")
res_json = res.json()
smiles = res_json["data"]["relationships"]["versions"]["data"][0]["attributes"]["smiles"]
mol = Chem.MolFromSmiles(smiles)
fp = AllChem.Get