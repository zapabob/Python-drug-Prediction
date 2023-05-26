from chembl_webresource_client.new_client import new_client
from keras.models import Sequential
from keras.layers import Dense
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import numpy as np

# ChEMBLデータベースからデータを取得
target = new_client.target
target_query = target.search('CHEMBL238')
targets = target_query.all()

activity = new_client.activity
res = activity.filter(target_chembl_id=targets[0]['target_chembl_id']).filter(standard_type="IC50")

# データの前処理
mols = [Chem.MolFromSmiles(compound['canonical_smiles']) for compound in res]
fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mols]
X = np.array(fps)
y = np.array([float(compound['standard_value']) for compound in res])

# 深層学習モデルの訓練
model = Sequential()
model.add(Dense(1024, input_dim=1024, activation='relu'))
model.add(Dense(512, activation='relu'))
model.add(Dense(1, activation='linear'))

model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
model.fit(X, y, epochs=10, batch_size=32)
