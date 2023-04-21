#!/bin/env python3

import cirpy

class SubstanceName:
    def __init__(self, smiles: str | None = None, iupac: str | None = None) -> None:
        if smiles is not None and iupac is not None:
            raise Exception("Both smiles and iupac cannot be provided.")
        elif smiles is not None:
            self.converted_name = cirpy.resolve(smiles, representation='iupac_name')
        elif iupac is not None:
            self.converted_name = cirpy.resolve(iupac, representation='smiles')
        else:
            raise Exception("Either smiles or iupac must be provided.")

    def get(self) -> str:
        return self.converted_name
