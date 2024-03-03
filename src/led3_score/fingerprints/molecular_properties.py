import numpy
from rdkit.Chem.QED import qed
from rdkit.Chem.GraphDescriptors import BertzCT
from rdkit.Chem import Descriptors as desc, Crippen, AllChem, Lipinski


class MolecularProperties:

    # reimplementation of the molecular properties from drug ex

    def __init__(self):
        self.impelemented_properties = {'MW': desc.MolWt,
                    'logP': Crippen.MolLogP,
                    'HBA': AllChem.CalcNumLipinskiHBA,
                    'HBD': AllChem.CalcNumLipinskiHBD,
                    'Rotable': AllChem.CalcNumRotatableBonds,
                    'Amide': AllChem.CalcNumAmideBonds,
                    'Bridge': AllChem.CalcNumBridgeheadAtoms,
                    'Hetero': AllChem.CalcNumHeteroatoms,
                    'Heavy': Lipinski.HeavyAtomCount,
                    'Spiro': AllChem.CalcNumSpiroAtoms,
                    'FCSP3': AllChem.CalcFractionCSP3,
                    'Ring': Lipinski.RingCount,
                    'Aliphatic': AllChem.CalcNumAliphaticRings,
                    'Aromatic': AllChem.CalcNumAromaticRings,
                    'Saturated': AllChem.CalcNumSaturatedRings,
                    'HeteroR': AllChem.CalcNumHeterocycles,
                    'TPSA': AllChem.CalcTPSA,
                    'Valence': desc.NumValenceElectrons,
                    'MR': Crippen.MolMR,
                    'QED': qed,
                    'Bertz': BertzCT}

    # calculate the molecular properties that are specified in a list
    def get_properties_list(self, molecule, properties_list):
        properties = {}
        for property_name in properties_list:
            properties[property_name] = self.impelemented_properties[property_name](molecule)
        return properties

    # calculate the molecular properties that are specified in a list and return only the values
    def get_properties_list_values(self, molecule, properties_list) -> numpy.array:
        properties = {}
        for property_name in properties_list:
            properties[property_name] = self.impelemented_properties[property_name](molecule)
        
        calculated_properties = [properties[property_name] for property_name in properties_list]

        #transform calculated_properties to a nummpy array
        calculated_properties = numpy.array(calculated_properties)

        return calculated_properties



