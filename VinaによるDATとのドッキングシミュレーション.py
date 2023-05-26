import subprocess

def run_docking(receptor_pdb, ligand_sdf, output_file):
    # AutoDock Vinaのパスを指定します。
    vina = "/path/to/vina"

    # ドッキングパラメータを指定します。
    center_x = 0.0  # レセプターの中心のX座標
    center_y = 0.0  # レセプターの中心のY座標
    center_z = 0.0  # レセプターの中心のZ座標
    size_x = 20  # ドッキングボックスのX方向のサイズ
    size_y = 20  # ドッキングボックスのY方向のサイズ
    size_z = 20  # ドッキングボックスのZ方向のサイズ

    # AutoDock Vinaを実行します。
    subprocess.call([
        vina,
        "--receptor", receptor_pdb,
        "--ligand", ligand_sdf,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--out", output_file
    ])

# レセプターとリガンドのファイルを指定します。
receptor_pdb = "4x2b.pdb"
ligand_sdf = "Structure2D_compound_CID_165361748.sdf"

# 出力ファイル名を指定します。
output_file = "docking_results.pdbqt"

# ドッキングシミュレーションを実行します。
run_docking(receptor_pdb, ligand_sdf, output_file)
