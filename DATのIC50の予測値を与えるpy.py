# Required Libraries
import tkinter as tk
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np

# Data
data = {
    "Test drug": ["Cocaine", "Amphetamine", "Mephedrone", "Methylone", "MDPV", "S-MDPV", "R-MDPV", "3,4-Catechol-PV", "4-OH-3-MeO-PV", "α-PVP", "α-PBP", "α-PPP"],
    "DAT uptake inhibition IC50 (nM)": [211, 93, 762, 1323, 4.1, 2.1, 382, 11, 784, 12, 63, 197],
    "DAT/SERT ratio": [1.48, 36.75, 0.55, 0.77, 816.82, None, None, 900, 12, 833, 159, 50]
}

df = pd.DataFrame(data)

# Cleaning the data - Removing None values
df = df.dropna()

# Split the Data
X = df["DAT/SERT ratio"].values.reshape(-1,1)
y = df["DAT uptake inhibition IC50 (nM)"].values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Model
model = LinearRegression()

# Fit the Model
model.fit(X_train, y_train)

# GUI
window = tk.Tk()
window.title('DAT IC50 Predictor')

# Entry Field
entry_field = tk.Entry(window)
entry_field.pack()

# Function to Predict DAT IC50
def predict_dat_ic50():
    iupac = entry_field.get()
    dat_sert_ratio = df.loc[df['Test drug'] == iupac, 'DAT/SERT ratio']
    if len(dat_sert_ratio) == 0:
        prediction = "IUPAC name not found in database"
    else:
        prediction = model.predict(np.array([dat_sert_ratio]).reshape(-1,1))
        prediction = str(prediction[0])
    result_label.config(text=prediction)

# Predict Button
predict_button = tk.Button(window, text='Predict DAT IC50', command=predict_dat_ic50)
predict_button.pack()

# Result Label
result_label = tk.Label(window, text='')
result_label.pack()

window.mainloop()
