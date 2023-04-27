#!/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, PDBWriter
from Bio.PDB import PDBParser, PDBIO, Select, Residue, PDBList
import os
from vina import Vina

class GeneratePDB:
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
class GetPDB:
    def __init__(self, pdb_id: str) -> None:
        self.pdb_id = pdb_id
        self.pdb_list = PDBList()
        self.pdb_list.retrieve_pdb_file(self.pdb_id, pdir=".", file_format="pdb", overwrite=True)

        self.remove_ligand_from_pdb(f"pdb{self.pdb_id.lower()}.ent")

    # remove ligands
    def remove_ligand_from_pdb(self, filename: str) -> None:
        mol = Chem.MolFromPDBFile(filename, removeHs=False, sanitize=False)
        new_mol = Chem.EditableMol(Chem.Mol())

        for atom in mol.GetAtoms():
            if atom.GetPDBResidueInfo().GetIsHeteroAtom() == False:
                new_mol.AddAtom(atom)

        final_mol = new_mol.GetMol()

        with PDBWriter(filename) as writer:
            writer.write(final_mol)


# class StardardResidues(Select):
#     def accept_residue(self, residue: Residue) -> bool:
#         return residue.get_resname() in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
#                                          "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
#                                          "TYR", "VAL"]

# # remove non-standard residues
# def remove_non_standard_residues(input_pdb: str, output_pdb: str) -> None:
#     parser = PDBParser()
#     structure = parser.get_structure("protein", input_pdb)
#     io = PDBIO()
#     io.set_structure(structure)
#     io.save(output_pdb, StardardResidues())


# generate the legand pdb file
def generate_ligand_pdb(smiles: str, filename: str, write_image: bool = False) -> None:
    # make the instance of the class
    generate_pdb = GeneratePDB(smiles)

    # draw the molecule
    if write_image:
        generate_pdb.draw(f"{filename}.png")

    # write the pdb file
    generate_pdb.write(f"{filename}.pdbqt")

    # post-process the pdb file
    # replace UNL with 3 spaces
    result: list[str] = []
    with open(f"{filename}.pdb", "r") as f:
        result = f.readlines()
    result = [line.replace("UNL", "   ") for line in result]
    with open(f"{filename}.pdb", "w") as f:
        f.writelines(result)


# genrate the receptor pdb file
def generate_receptor_pdb(pdb_id: str, filename: str) -> None:
    # get the pdb file
    GetPDB(pdb_id) # Dopamine Transporter

    # rename the pdb file
    os.rename(f"pdb{pdb_id.lower()}.ent", f"{filename}.pdbqt")


# the main function to test the class
def main() -> None:
    LIGAND_FILE_NAME: str = "02_meth"
    LIGAND_SMILES: str = "CNC(C)Cc1ccccc1"

    RECEPTOR_FILE_NAME: str = "02_Receptor"
    RECEPTOR_PDB_ID: str = "4M48"

    # generate the regand pdbqt file
    generate_ligand_pdb(LIGAND_SMILES, LIGAND_FILE_NAME)

    # generate the receptor pdb file
    generate_receptor_pdb(RECEPTOR_PDB_ID, RECEPTOR_FILE_NAME)

    # docking with vina
    print("Simulation Start")
    v = Vina(sf_name="vina", cpu=16)
    v.set_receptor(f"{RECEPTOR_FILE_NAME}.pdbqt")
    v.set_ligand_from_file(f"{LIGAND_FILE_NAME}.pdbqt")
    v.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
    v.dock(exhaustiveness=32, out="out.pdbqt")
    v.write_poses("poses.pdbqt", n_poses=5, overwrite=True)
    print("Simulation End")


if __name__ == "__main__":
    main()
