#!/bin/env python3

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, SDWriter


class GeneratePDB:
    def __init__(self, smiles: str) -> None:
        self.mol = Chem.MolFromSmiles(smiles)
        Chem.AddHs(self.mol)
        AllChem.Compute2DCoords(self.mol)

    def draw(self, filename: str) -> None:
        Draw.MolsToGridImage([self.mol], molsPerRow=1, subImgSize=(500, 500)).save(filename)

    def write(self, file_name: str) -> None:
        Chem.rdmolfiles.MolToPDBFile(self.mol, file_name)


# the main function to test the class
def main() -> None:
    # make the instance of the class
    generate_pdb = GeneratePDB("CNC(C)Cc1ccccc1")

    # draw the molecule
    # generate_pdb.draw("02_meth.png")

    # write the pdb file
    generate_pdb.write("02_meth.pdb")


if __name__ == "__main__":
    main()
