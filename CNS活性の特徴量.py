import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

# Step 1: データの準備
data = pd.read_csv('MA.csv','MPH.csv','Cocaine.csv')  # "your_data.csv"は実際のファイルパスに置き換えてください

# 欠損値処理（例として平均値で補完）
data = data.fillna(data.mean())

# Step 2: 特徴量の抽出
# この例では、すべてのカラムが数値データで、"CNS_activity"が目的変数で、それ以外が特徴量であると仮定しています。
X = data.drop("CNS_activity", axis=1)
y = data["CNS_activity"]

# Step 3: モデルの訓練
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
model = LinearRegression()
model.fit(X_train, y_train)

# Step 4: モデルの評価
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
print(f"Mean Squared Error: {mse}")
