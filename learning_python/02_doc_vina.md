# AutoDock Vinaの使用方法

## SMILES記法からPDBファイルを生成
```bash
./02_vina.py
```
これで**02_meth.pdb**というファイルができる。

## AutoDock Vinaのセットアップ
```bash
AutoDock Vina3_setup.py -s 400 -g 200 --noxt --noh 02_Receptor.pdb 02_meth.pdb
```

## ドッキングシミュレーションの実行
```bash
aaa
```

## 中間ファイルの消し方
```bash
./02_clean.sh
```