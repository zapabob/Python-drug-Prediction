import pandas as pd
import numpy as np
import tkinter as tk
from rdkit import Chem
import deepchem as dc
from deepchem.models import GraphConvModel
import requests
import time

# Prepare and preprocess data
compounds = ["MDMA", "Methylone", "BDB", "MBDB", "MDPV", "2C-I", "2C-E", "2C-C", "TMA", "TMA-2", "TMA-6", "PMMA", "4FMP", "MAP"]
data = {
    "Compound": compounds,
    "DA IC50": [1.4e-6, 2.9e-6, 7.9e-6, 6.3e-6, 4.3e-8, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 1.4e-5, 7.7e-7, 3.7e-7],
}

df = pd.DataFrame(data)
df['DA IC50'] = np.log(df['DA IC50'])

# Get SMILES representations from PubChem
def get_smiles(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
    else:
        return np.nan

df['SMILES'] = [get_smiles(compound) for compound in compounds]
df = df.dropna(subset=['SMILES'])  # Drop rows with no SMILES

# Define the featurizer
featurizer = dc.feat.ConvMolFeaturizer()

# Load the dataset
loader = dc.data.CSVLoader(tasks=['DA IC50'], smiles_field="SMILES", featurizer=featurizer)
dataset = loader.featurize(df)

# Split the dataset
splitter = dc.splits.RandomSplitter()
train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(dataset)

# Train the model
model = GraphConvModel(n_tasks=1, mode='regression', dropout=0.2)
model.fit(train_dataset, nb_epoch=100)

def predict_ic50(smiles):
    dataset = loader.featurize(pd.DataFrame({"SMILES": [smiles]}))
    prediction = model.predict(dataset)
    return np.exp(prediction[0][0])

def predict():
    try:
        compound = compound_entry.get()
        smiles = get_smiles(compound)
        if pd.isnull(smiles):
            ic50_text.set(f"No SMILES found for {compound}")
        else:
                        ic50 = predict_ic50(smiles)
        ic50_text.set(f"Predicted IC50: {ic50}")
    except Exception as e:
        ic50_text.set(f"Error: {str(e)}")

root = tk.Tk()
root.title("IC50 Predictor")
root.geometry("300x200")

ic50_text = tk.StringVar()

compound_label = tk.Label(root, text="Enter Compound:")
compound_label.pack()
compound_entry = tk.Entry(root)
compound_entry.pack()

predict_button = tk.Button(root, text="Predict IC50", command=predict)
predict_button.pack()

ic50_label = tk.Label(root, textvariable=ic50_text)
ic50_label.pack()

root.mainloop()
