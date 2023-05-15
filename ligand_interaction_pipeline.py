#!/bin/env python3
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
# データの読み込み
data = pd.read_csv('A.csv')

# 前処理（欠損値補完、外れ値除去、カテゴリ変数のエンコーディングなど）
# ...
#欠損値の補完
data.fillna(data.mean(), inplace=True)
#カテゴリ変数のエンコーディング
data = pd.get_dummies(data, drop_first=True)
#外れ値の処理
#外れ値の処理（IQR法を用いた外れ値の検出と削除）
Q1 = data.quantile(0.25)
Q3 = data.quantile(0.75)
IQR = Q3 - Q1

data = data[~((data < (Q1 - 1.5 * IQR)) | (data > (Q3 + 1.5 * IQR))).any(axis=1)]

# データの分割
X = data.drop('target', axis=1)
y = data['target']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 特徴量スケーリング
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# モデルの選択（例: ランダムフォレスト）
model = RandomForestClassifier()

# ハイパーパラメータチューニング
param_grid = {
'n_estimators': [10, 50, 100, 200],
'max_depth': [None, 10, 20, 30],
'min_samples_split': [2, 5, 10]
}
grid_search = GridSearchCV(model, param_grid, cv=5)
grid_search.fit(X_train_scaled, y_train)

# 最適なハイパーパラメータでモデルを再学習
best_model = grid_search.best_estimator_
best_model.fit(X_train_scaled, y_train)

# テストデータでの評価
y_pred = best_model.predict(X_test_scaled)
accuracy = accuracy_score(y_test, y_pred)
print(f'Test accuracy: {accuracy}')
