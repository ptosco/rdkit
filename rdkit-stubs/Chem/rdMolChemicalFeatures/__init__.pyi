"""Module containing from chemical feature and functions to generate the"""
from __future__ import annotations
import rdkit.Chem.rdMolChemicalFeatures
import typing
import Boost.Python

__all__ = [
    "BuildFeatureFactory",
    "BuildFeatureFactoryFromString",
    "GetAtomMatch",
    "MolChemicalFeature",
    "MolChemicalFeatureFactory"
]


class MolChemicalFeature(Boost.Python.instance):
    """
    Class to represent a chemical feature.
        These chemical features may or may not have been derived from molecule object;
        i.e. it is possible to have a chemical feature that was created just from its type
        and location.
    """
    @staticmethod
    def ClearCache( arg1: MolChemicalFeature) -> None: 
        """
        ClearCache( arg1: MolChemicalFeature) -> None
            Clears the cache used to store position information.

            C++ signature :
                void ClearCache(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def GetActiveConformer( arg1: MolChemicalFeature) -> int: 
        """
        GetActiveConformer( arg1: MolChemicalFeature) -> int
            Gets the conformer to use.

            C++ signature :
                int GetActiveConformer(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def GetAtomIds( arg1: MolChemicalFeature) -> object: 
        """
        GetAtomIds( arg1: MolChemicalFeature) -> object
            Get the IDs of the atoms that participate in the feature

            C++ signature :
                _object* GetAtomIds(RDKit::MolChemicalFeature)
        """
    @staticmethod
    def GetFactory( arg1: MolChemicalFeature) -> MolChemicalFeatureFactory: 
        """
        GetFactory( arg1: MolChemicalFeature) -> MolChemicalFeatureFactory
            Get the factory used to generate this feature

            C++ signature :
                RDKit::MolChemicalFeatureFactory const* GetFactory(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def GetFamily( arg1: MolChemicalFeature) -> str: 
        """
        GetFamily( arg1: MolChemicalFeature) -> str
            Get the family to which the feature belongs; donor, acceptor, etc.

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetFamily(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def GetId( arg1: MolChemicalFeature) -> int: 
        """
        GetId( arg1: MolChemicalFeature) -> int
            Returns the identifier of the feature
            

            C++ signature :
                int GetId(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def GetMol( arg1: MolChemicalFeature) -> Mol: 
        """
        GetMol( arg1: MolChemicalFeature) -> Mol
            Get the molecule used to derive the features

            C++ signature :
                RDKit::ROMol const* GetMol(RDKit::MolChemicalFeature {lvalue})
        """
    @typing.overload
    def GetPos(self, confId: int) -> Point3D: 
        """
        GetPos( self: MolChemicalFeature, confId: int) -> Point3D
            Get the location of the chemical feature

            C++ signature :
                RDGeom::Point3D GetPos(RDKit::MolChemicalFeature {lvalue},int)

            C++ signature :
                RDGeom::Point3D GetPos(RDKit::MolChemicalFeature {lvalue})
        """
    @typing.overload
    def GetPos(self) -> Point3D: ...
    @staticmethod
    def GetType( arg1: MolChemicalFeature) -> str: 
        """
        GetType( arg1: MolChemicalFeature) -> str
            Get the specific type for the feature

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetType(RDKit::MolChemicalFeature {lvalue})
        """
    @staticmethod
    def SetActiveConformer( arg1: MolChemicalFeature, arg2: int) -> None: 
        """
        SetActiveConformer( arg1: MolChemicalFeature, arg2: int) -> None
            Sets the conformer to use (must be associated with a molecule).

            C++ signature :
                void SetActiveConformer(RDKit::MolChemicalFeature {lvalue},int)
        """
    pass
class MolChemicalFeatureFactory(Boost.Python.instance):
    """
    Class to featurize a molecule
    """
    @staticmethod
    def GetFeatureDefs( arg1: MolChemicalFeatureFactory) -> dict: 
        """
        GetFeatureDefs( arg1: MolChemicalFeatureFactory) -> dict
            Get a dictionary with SMARTS definitions for each feature type

            C++ signature :
                boost::python::dict GetFeatureDefs(RDKit::MolChemicalFeatureFactory)
        """
    @staticmethod
    def GetFeatureFamilies( arg1: MolChemicalFeatureFactory) -> tuple: 
        """
        GetFeatureFamilies( arg1: MolChemicalFeatureFactory) -> tuple
            Get a tuple of feature types

            C++ signature :
                boost::python::tuple GetFeatureFamilies(RDKit::MolChemicalFeatureFactory)
        """
    @staticmethod
    def GetMolFeature( arg1: MolChemicalFeatureFactory, mol: Mol, idx: int, includeOnly: str = '', recompute: bool = True, confId: int = -1) -> MolChemicalFeature: 
        """
        GetMolFeature( arg1: MolChemicalFeatureFactory, mol: Mol, idx: int, includeOnly: str = '', recompute: bool = True, confId: int = -1) -> MolChemicalFeature
            returns a particular feature (by index)

            C++ signature :
                boost::shared_ptr<RDKit::MolChemicalFeature> GetMolFeature(RDKit::MolChemicalFeatureFactory,RDKit::ROMol,int [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,bool=True [,int=-1]]])
        """
    @staticmethod
    def GetNumFeatureDefs( arg1: MolChemicalFeatureFactory) -> int: 
        """
        GetNumFeatureDefs( arg1: MolChemicalFeatureFactory) -> int
            Get the number of feature definitions

            C++ signature :
                int GetNumFeatureDefs(RDKit::MolChemicalFeatureFactory {lvalue})
        """
    @staticmethod
    def GetNumMolFeatures( arg1: MolChemicalFeatureFactory, mol: Mol, includeOnly: str = '') -> int: 
        """
        GetNumMolFeatures( arg1: MolChemicalFeatureFactory, mol: Mol, includeOnly: str = '') -> int
            Get the number of features the molecule has

            C++ signature :
                int GetNumMolFeatures(RDKit::MolChemicalFeatureFactory,RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
        """
    pass
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
