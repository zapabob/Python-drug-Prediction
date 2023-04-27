#!/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, PDBWriter
from Bio.PDB import PDBParser, PDBIO, Select, Residue, PDBList
import os
from vina import Vina


class Ligand:
    def __init__(self, smiles: str) -> None:
        self.mol = Chem.MolFromSmiles(smiles)
        self.mol = Chem.AddHs(self.mol)
        p = AllChem.ETKDGv2()
        AllChem.EmbedMolecule(self.mol, p)


    def draw(self, filename: str) -> None:
        Draw.MolsToGridImage([self.mol], molsPerRow=1, subImgSize=(500, 500)).save(filename)


    def write(self, file_name: str) -> None:
        Chem.rdmolfiles.MolToPDBFile(self.mol, file_name)


# The class to get PDB file from PDB ID
class Receptor:
    def __init__(self, pdb_id: str) -> None:
        self.pdb_id = pdb_id
        self.pdb_list = PDBList()


    # write the cif file
    def write(self, filename: str) -> None:
        self.pdb_list.retrieve_pdb_file(self.pdb_id, pdir=".", file_format="mmCif", overwrite=True)

        # rename the file
        os.rename(f"{self.pdb_id.lower()}.cif", f"{filename}.cif")


# generate the legand pdb file
def generate_ligand_pdb(smiles: str, filename: str, write_image: bool = False) -> None:
    # make the instance of the class
    lig = Ligand(smiles)

    # draw the molecule
    if write_image:
        lig.draw(f"{filename}.png")

    # write the pdb file
    lig.write(f"{filename}.pdb")

    # post-process the pdb file
    # replace UNL with 3 spaces
    result: list[str] = []
    with open(f"{filename}.pdb", "r") as f:
        result = f.readlines()
    result = [line.replace("UNL", "   ") for line in result]
    with open(f"{filename}.pdb", "w") as f:
        f.writelines(result)


# genrate the receptor cif file
def generate_receptor_cif(pdb_id: str, filename: str) -> None:
    # get the pdb file
    Receptor(pdb_id).write(filename)


# the main function to test the class
def main() -> None:
    LIGAND_FILE_NAME: str = "02_meth"
    LIGAND_SMILES: str = "CNC(C)Cc1ccccc1"

    RECEPTOR_FILE_NAME: str = "02_Receptor"
    RECEPTOR_PDB_ID: str = "4XP1"

    # generate the regand pdb file
    generate_ligand_pdb(LIGAND_SMILES, LIGAND_FILE_NAME)

    # generate the receptor pdb file
    generate_receptor_cif(RECEPTOR_PDB_ID, RECEPTOR_FILE_NAME)

    # # docking with vina
    # print("Simulation Start")
    # v = Vina(sf_name="vina", cpu=16)
    # v.set_receptor(f"{RECEPTOR_FILE_NAME}.pdbqt")
    # v.set_ligand_from_file(f"{LIGAND_FILE_NAME}.pdbqt")
    # v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
    # v.dock(exhaustiveness=32, out="out.pdbqt")
    # v.write_poses("poses.pdbqt", n_poses=5, overwrite=True)
    # print("Simulation End")


if __name__ == "__main__":
    main()
