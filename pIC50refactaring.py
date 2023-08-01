import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score, roc_curve
import statsmodels.api as sm
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

# Extracting data from the provided file
with open(.py) as file:
    exec(file.read())
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
    "pemoline":'C1=CC=C(C=C1)C2C(=O)N=C(O2)N',
    "dopamine":'NCCc1ccc(O)c(O)c1',
}
# Create a DataFrame with estimated pIC50 values
data = pd.DataFrame(list(compounds.values()), columns=["SMILES", "pIC50"], index=compounds.keys())

# Compute LogP and LogS for each compound
data["LogP"] = data["SMILES"].apply(lambda x: Descriptors.MolLogP(Chem.MolFromSmiles(x)))
data["LogS"] = data["SMILES"].apply(lambda x: Descriptors.MolLogS(Chem.MolFromSmiles(x)))

# Define X (predictors) and y (response)
X = data[["LogP", "LogS"]]
y = data["pIC50"]

# Fit logistic regression model
log_reg = LogisticRegression(max_iter=1000)
log_reg.fit(X, y)

# Predict the probabilities
data["PredictedProb"] = log_reg.predict_probate(X)[:,1]

# Model Summary using statsmodels
logits_model = sm.Logits(y, sm.add_constant(X))
result = logits_model.fit()
print(result.summary())

# ROC Curve
fpr, tpr, thresholds = roc_curve(y, data["PredictedProb"])
plt.plot(fpr, tpr, label=f'AUC: {roc_auc_score(y, data["PredictedProb"]):.2f}')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend(loc='lower right')
plt.show()

