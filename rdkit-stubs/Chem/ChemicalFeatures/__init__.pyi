from __future__ import annotations
import rdkit.Chem.ChemicalFeatures
import typing
from rdkit.Chem.rdChemicalFeatures import FreeChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeatureFactory

__all__ = [
    "BuildFeatureFactory",
    "BuildFeatureFactoryFromString",
    "FreeChemicalFeature",
    "GetAtomMatch",
    "MCFF_GetFeaturesForMol",
    "MolChemicalFeature",
    "MolChemicalFeatureFactory"
]


def BuildFeatureFactory( arg1: str) -> MolChemicalFeatureFactory:
    """
    BuildFeatureFactory( arg1: str) -> MolChemicalFeatureFactory
        Construct a feature factory given a feature definition in a file

        C++ signature :
            RDKit::MolChemicalFeatureFactory* BuildFeatureFactory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def BuildFeatureFactoryFromString( arg1: str) -> MolChemicalFeatureFactory:
    """
    BuildFeatureFactoryFromString( arg1: str) -> MolChemicalFeatureFactory
        Construct a feature factory given a feature definition block

        C++ signature :
            RDKit::MolChemicalFeatureFactory* BuildFeatureFactoryFromString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def GetAtomMatch( featMatch: AtomPairsParameters, maxAts: int = 1024) -> object:
    """
    GetAtomMatch( featMatch: AtomPairsParameters, maxAts: int = 1024) -> object
        Returns an empty list if any of the features passed in share an atom.
         Otherwise a list of lists of atom indices is returned.
        

        C++ signature :
            boost::python::api::object GetAtomMatch(boost::python::api::object [,int=1024])
    """
