import cirpy

def iupac_to_smiles(iupac):
    return cirpy.resolve(iupac, "smiles")

A=iupac_to_smiles("1-methyl-phenylpropaneamine")
print(A)
