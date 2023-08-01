from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

class TGDDistanceCalculator:
    PHENETHYLAMINE_SMILES = "C1=CC=C(C=C1)CCCN"
    
    def __init__(self, compounds):
        self.compounds = compounds
        self.reference_distance = self._compute_reference_distance()

    def _identify_atoms(self, mol):
        aromatic_carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic() and a.GetSymbol() == 'C']
        amine_nitrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'N' and a.GetTotalDegree() == 3]
        return aromatic_carbons, amine_nitrogens

    def _generate_3D_coordinates(self, mol):
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        return mol

    def _compute_reference_distance(self):
        mol = Chem.MolFromSmiles(self.PHENETHYLAMINE_SMILES)
        mol = self._generate_3D_coordinates(mol)
        aromatic_carbons, amine_nitrogens = self._identify_atoms(mol)
        return np.linalg.norm(
            np.array(mol.GetConformer().GetAtomPosition(aromatic_carbons[0])) -
            np.array(mol.GetConformer().GetAtomPosition(amine_nitrogens[0]))
        )

    def compute_normalized_distance(self, mol):
        mol = self._generate_3D_coordinates(mol)
        aromatic_carbons, amine_nitrogens = self._identify_atoms(mol)
        distance = np.linalg.norm(
            np.array(mol.GetConformer().GetAtomPosition(aromatic_carbons[0])) -
            np.array(mol.GetConformer().GetAtomPosition(amine_nitrogens[0]))
        )
        return 100 * (distance / self.reference_distance)

    def compute_all_normalized_distances(self):
        distances = {}
        for name, smiles in self.compounds.items():
            mol = Chem.MolFromSmiles(smiles)
            distances[name] = self.compute_normalized_distance(mol)
        return distances

# Define compounds
compounds = {...}  # Your compounds dictionary here

# Compute and display normalized distances
calculator = TGDDistanceCalculator(compounds)
distances = calculator.compute_all_normalized_distances()
for name, distance in distances.items():
    print(f"Normalized distance for {name}: {distance}")
