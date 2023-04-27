#!/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from Bio.PDB import PDBParser, PDBIO, Select, Residue, PDBList
import os


class GeneratePDB:
    def __init__(self, smiles: str) -> None:
        self.mol = Chem.MolFromSmiles(smiles)
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


class StardardResidues(Select):
    def accept_residue(self, residue: Residue) -> bool:
        return residue.get_resname() in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                                         "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                                         "TYR", "VAL"]

# remove non-standard residues
def remove_non_standard_residues(input_pdb: str, output_pdb: str) -> None:
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, StardardResidues())


# the main function to test the class
def main() -> None:
    # make the instance of the class
    generate_pdb = GeneratePDB("CNC(C)Cc1ccccc1")

    # draw the molecule
    # generate_pdb.draw("02_meth.png")

    # write the pdb file
    generate_pdb.write("02_meth.pdb")

    # get the pdb file
    PDB_ID: str = "4M48"
    GetPDB(PDB_ID) # Dopamine Transporter

    # remove non-standard residues
    filename: str = f"pdb{PDB_ID.lower()}.ent"
    remove_non_standard_residues(filename, "02_Receptor.pdb")

    # remove pdb4m48.ent
    os.remove(filename)


if __name__ == "__main__":
    main()
