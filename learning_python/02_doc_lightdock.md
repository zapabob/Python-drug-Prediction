# LightDockの使用方法

## SMILES記法からPDBファイルを生成
```bash
./02_lightdock.py
```
これで**02_meth.pdb**というファイルができる。

## LightDockのセットアップ
```bash
lightdock3_setup.py -s 400 -g 200 --noxt --noh 02_DAT.pdb 02_meth.pdb
```

## ドッキングシミュレーションの実行
```bash
aaa
```

## 中間ファイルの消し方
```bash
rm -rf init swarm* setup.json lightdock*
```