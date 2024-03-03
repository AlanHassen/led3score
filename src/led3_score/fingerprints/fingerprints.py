from abc import ABC, abstractmethod
from led3_score.fingerprints.molecular_properties import MolecularProperties

import numpy
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from loguru import logger



class Fingerprint(ABC):
    """
    Abstract class for fingerprints.
    """

    def __init__(self, fingerprint_name:str) -> None:
        self.fingerprint_name = fingerprint_name
        logger.info(f"Fingerprint {fingerprint_name} initialized.")

    @abstractmethod
    def create_fingerprint(smiles:str):
        """
        Creates the fingerprint.
        """
        raise NotImplementedError


class RDKitBasedFingerprint(Fingerprint):

    """Class for creating fingerprints using RDKit"""

    @staticmethod
    def smiles_to_mol(smiles):
        """
        Transform a smiles string into a Rdkit mol
        Args:
            smiles: the smiles string

        Returns:
            rdkit mol
        """
        molecule = Chem.MolFromSmiles(smiles)
        return molecule

    @staticmethod
    def numpy_to_rdkit_datastructure(numpy_array):
        """
        Transform the numpy array of a fingerprint into a Rdkit datastructure

        Args:
            numpyArray: Numpy array representing a molecular fingerprint

        Returns:
            RDkit fingerprint datastructure
        """
        # Note: we need to replace the trues with ints and than transform that into strings
        int_array = numpy_array.astype(int)
        bit_string = "".join(int_array.astype(str))
        rdkit_fingerprint = DataStructs.cDataStructs.CreateFromBitString(bit_string)
        return rdkit_fingerprint

    @staticmethod
    def rdkit_datastructure_to_numpy(rdkit_datastructure):
        """
        Transform a RDkit fingerprint datastructure to a numpy array

        Returns:
            numpy_array of the fingerprint
        """
        numpy_array = numpy.zeros((0,), dtype=numpy.bool_)
        DataStructs.ConvertToNumpyArray(rdkit_datastructure, numpy_array)
        return numpy_array

    def __init__(self, fingerprint_name:str, transform_to_numpy_array: bool) -> None:
        self.transform_to_numpy_array = transform_to_numpy_array
        super().__init__(fingerprint_name = fingerprint_name)

    def create_fingerprint(self, smiles):
        """
        Create a binary fingerprint from a smiles
        Args:
            smiles: the molecule smiles

        Returns:
            binary fingerprint
        """

        molecule = self.smiles_to_mol(smiles)
        fingerprint = self.create_rdkit_based_fingerprint(molecule)
        if self.transform_to_numpy_array:
            return self.rdkit_datastructure_to_numpy(fingerprint)
        else:
            return fingerprint

    @abstractmethod
    def create_rdkit_based_fingerprint(self, molecule):
        raise NotImplementedError


class MorganFingerprint(RDKitBasedFingerprint):
    """
    Abstract class for fingerprints based on Morgan fingerprints.
    """
    def __init__(self, fingerprint_name: str, transform_to_numpy_array: bool, radius:int, number_of_bits:int, use_features: bool) -> None:
        self.radius = radius
        self.number_of_bits = number_of_bits
        self.use_features = use_features
        super().__init__(fingerprint_name, transform_to_numpy_array)

    @abstractmethod
    def create_rdkit_based_fingerprint(self, molecule):
        raise NotImplementedError

class HashedMorganFingerprint(MorganFingerprint):
    """
    Fingerprint based on Hashed Morgan fingerprints.
    """

    def __init__(self, transform_to_numpy_array: bool, radius:int, number_of_bits:int, use_features: bool) -> None:
        fingerprint_name = f"Hashed Morgan Fingerprint with Radius {radius} and {number_of_bits} bits using Features: {use_features}"
        super().__init__(fingerprint_name, transform_to_numpy_array, radius, number_of_bits, use_features)

    def create_rdkit_based_fingerprint(self, molecule):
        ecfp_hashed_fingerprint = AllChem.GetHashedMorganFingerprint(molecule, radius=self.radius, nBits=self.number_of_bits, useFeatures=self.use_features)
        return ecfp_hashed_fingerprint


class BinaryMorganFingerprint(MorganFingerprint):

    """
    Fingerprint based on Binary Morgan fingerprints.
    """

    def __init__(self, transform_to_numpy_array: bool, radius:int, number_of_bits:int, use_features: bool) -> None:
        fingerprint_name = f"Binary Morgan Fingerprint with Radius {radius} and {number_of_bits} bits using Features: {use_features}"
        super().__init__(fingerprint_name, transform_to_numpy_array, radius, number_of_bits, use_features)

    def create_rdkit_based_fingerprint(self, molecule):
        ecfp_binary_fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=self.radius, nBits=self.number_of_bits, useFeatures=self.use_features)
        return ecfp_binary_fingerprint

class CountBasedMorganFingerprint(MorganFingerprint):

    """
    Fingerprint based on Count Based Morgan fingerprints.
    """

    def __init__(self, transform_to_numpy_array: bool, radius:int, number_of_bits:int, use_features: bool) -> None:
        fingerprint_name = f"Count Based Morgan Fingerprint with Radius {radius} and {number_of_bits} bits using Features: {use_features}"
        assert transform_to_numpy_array == False, "Count Based Morgan Fingerprints don't need to be converted to numpy arrays. Therefore transform_to_numpy_array should be False."
        super().__init__(fingerprint_name, transform_to_numpy_array, radius, number_of_bits, use_features)

    def create_rdkit_based_fingerprint(self, molecule):
        """
        Implementation of the count based fingerprint with radius 3 and 2048 bits and features
        Source: https://github.com/reymond-group/RAscore/blob/master/RAscore/RAscore_XGB.py
        """
        ecfp_fingerprint = AllChem.GetMorganFingerprint(molecule, radius=self.radius, useFeatures=self.use_features)

        ecfp_count_fingerprint = numpy.zeros((self.number_of_bits,), numpy.int32)
        for idx, v in ecfp_fingerprint.GetNonzeroElements().items():
            nidx = idx % self.number_of_bits
            ecfp_count_fingerprint[nidx] += int(v)
        return ecfp_count_fingerprint

class DrugExFingerprint(MorganFingerprint):

    """Hardcoded version of the DrugEx fingerprint.
    Source: https://gitlab.services.universiteitleiden.nl/cdd/DrugEx/-/blob/dev/drugex/training/scorers/predictors.py
    """

    def __init__(self) -> None:
        transform_to_numpy_array = False
        radius = 3
        number_of_bits = 2048
        use_features = False
        fingerprint_name = f"DrugEx Fingerprint (Binary Morgan Fingerprint with Radius {radius} and {number_of_bits} bits using Features: {use_features} + Chemical Features)"
        assert transform_to_numpy_array == False, "Count Based Morgan Fingerprints don't need to be converted to numpy arrays. Therefore transform_to_numpy_array should be False."
        super().__init__(fingerprint_name, transform_to_numpy_array, radius, number_of_bits, use_features)

    def create_rdkit_based_fingerprint(self, molecule):
        ecfp_binary_fingerprint = AllChem.GetMorganFingerprintAsBitVect(molecule, radius=self.radius, nBits=self.number_of_bits, useFeatures=self.use_features)
        ecfp_numpy_array = self.rdkit_datastructure_to_numpy(ecfp_binary_fingerprint)
        chemical_properties_numpy_array = self._caluclate_chemical_properties(molecule)
        combined_fingerprint = numpy.concatenate([ecfp_numpy_array, chemical_properties_numpy_array])
        return combined_fingerprint

    def _caluclate_chemical_properties(self, molecule):
        """Chemical Properties added to the fingerprint in DrugEx
        Source: https://gitlab.services.universiteitleiden.nl/cdd/DrugEx/-/blob/dev/drugex/training/scorers/predictors.py

        """

        prop_list = ['MW', 'logP', 'HBA', 'HBD', 'Rotable', 'Amide',
                     'Bridge', 'Hetero', 'Heavy', 'Spiro', 'FCSP3', 'Ring',
                     'Aliphatic', 'Aromatic', 'Saturated', 'HeteroR', 'TPSA', 'Valence', 'MR']
        
        molecularProperties = MolecularProperties()
        
        chemical_descriptor_values = molecularProperties.get_properties_list_values(molecule, prop_list)

        return chemical_descriptor_values
    


