import tkinter as tk
from keras.models import load_model
from rdkit.Chem import MolFromSmiles, MolToSmiles
import numpy as np

# モデルの読み込み
models = {target_id: load_model(f'model_{target_id}.h5') for target_id in targets_ids}

# IC50を予測する関数
def predict_ic50(target_id, iupac_name):
    smiles = MolToSmiles(MolFromSmiles(iupac_name))
    descriptors = compute_descriptors(smiles)
    predicted_ic50 = models[target_id].predict(np.array([descriptors]))
    return predicted_ic50

# GUIの作成
root = tk.Tk()
root.title("Monoamine transporter IC50 Ratio Predictor")

# 入力フィールド
iupac_entry = tk.Entry(root)
iupac_entry.pack()

# 結果表示テキストボックス
result_textbox = tk.Text(root)
result_textbox.pack()

# ボタンが押されたときの処理
def on_predict_button_press():
    iupac_name = iupac_entry.get()
    try:
        ic50_dat = predict_ic50('CHEMBL240', iupac_name)
        ic50_sert = predict_ic50('CHEMBL228', iupac_name)
        ic50_net = predict_ic50('CHEMBL222', iupac_name)

        pic50_dat = -np.log10(ic50_dat)
        pic50_sert = -np.log10(ic50_sert)
        pic50_net = -np.log10(ic50_net)

        ratio_dat_sert = pic50_dat / pic50_sert
        ratio_dat_net = pic50_dat / pic50_net

        result_textbox.insert(tk.END, f"Dopamine transporter / Serotonin transporter pIC50 ratio: {ratio_dat_sert}\n"
                                    f"Dopamine transporter / Norepinephrine transporter pIC50 ratio: {ratio_dat_net}")
    except Exception as e:
        result_textbox.insert(tk.END, f"Error: {str(e)}")

# リセットボタンが押されたときの処理
def on_reset_button_press():
    iupac_entry.delete(0, tk.END)
    result_textbox.delete('1.0', tk.END)

# ボタン
predict_button = tk.Button(root, text="Predict pIC50 ratios", command=on_predict_button_press)
predict_button.pack()

# リセットボタン
reset_button = tk.Button(root, text="Reset", command=on_reset_button_press)
reset_button.pack()

root.mainloop()
