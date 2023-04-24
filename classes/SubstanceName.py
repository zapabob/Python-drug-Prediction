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
        # if self.converted_name is None, this method returns an empty string.
        return self.converted_name if self.converted_name is not None else ""


# the main function to test the class
def main() -> None:
    # make the instance of the class
    substance_name = SubstanceName(iupac='1,2,3,4-teoquinoline')

    # get the converted name
    converted_name = substance_name.get()

    # print the converted name
    print(converted_name)


if __name__ == "__main__":
    main()