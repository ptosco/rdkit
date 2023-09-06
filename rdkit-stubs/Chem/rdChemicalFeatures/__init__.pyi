"""Module containing free chemical feature functionality
     These are feature that are not associated with molecules. They are 
     are typically derived from pharmacophores and site-maps.
"""
from __future__ import annotations
import rdkit.Chem.rdChemicalFeatures
import typing
import Boost.Python

__all__ = [
    "FreeChemicalFeature"
]


class FreeChemicalFeature(Boost.Python.instance):
    """
    Class to represent a free chemical features.
        These chemical features are not associated with a molecule, though they can be matched 
        to molecular featufres
    """
    @staticmethod
    def GetFamily( arg1: FreeChemicalFeature) -> str: 
        """
        GetFamily( arg1: FreeChemicalFeature) -> str
            Get the family of the feature

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetFamily(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
    @staticmethod
    def GetId( arg1: FreeChemicalFeature) -> int: 
        """
        GetId( arg1: FreeChemicalFeature) -> int
            Get the id of the feature

            C++ signature :
                int GetId(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
    @staticmethod
    def GetPos( arg1: FreeChemicalFeature) -> Point3D: 
        """
        GetPos( arg1: FreeChemicalFeature) -> Point3D
            Get the position of the feature

            C++ signature :
                RDGeom::Point3D GetPos(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
    @staticmethod
    def GetType( arg1: FreeChemicalFeature) -> str: 
        """
        GetType( arg1: FreeChemicalFeature) -> str
            Get the sepcific type for the feature

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetType(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
    @staticmethod
    def SetFamily( arg1: FreeChemicalFeature, arg2: str) -> None: 
        """
        SetFamily( arg1: FreeChemicalFeature, arg2: str) -> None
            Set the family of the feature

            C++ signature :
                void SetFamily(ChemicalFeatures::FreeChemicalFeature {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetId( arg1: FreeChemicalFeature, arg2: int) -> None: 
        """
        SetId( arg1: FreeChemicalFeature, arg2: int) -> None
            Set the id of the feature

            C++ signature :
                void SetId(ChemicalFeatures::FreeChemicalFeature {lvalue},int)
        """
    @staticmethod
    def SetPos( arg1: FreeChemicalFeature, arg2: Point3D) -> None: 
        """
        SetPos( arg1: FreeChemicalFeature, arg2: Point3D) -> None
            Set the feature position

            C++ signature :
                void SetPos(ChemicalFeatures::FreeChemicalFeature {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def SetType( arg1: FreeChemicalFeature, arg2: str) -> None: 
        """
        SetType( arg1: FreeChemicalFeature, arg2: str) -> None
            Set the sepcific type for the feature

            C++ signature :
                void SetType(ChemicalFeatures::FreeChemicalFeature {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __getinitargs__( arg1: FreeChemicalFeature) -> tuple: 
        """
        __getinitargs__( arg1: FreeChemicalFeature) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(ChemicalFeatures::FreeChemicalFeature)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: 
        """
        __init__( arg1: object, arg2: str) -> None

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D [,int=-1])

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, family: str, type: str, loc: Point3D, id: int = -1) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, family: str, loc: Point3D) -> None: ...
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 136
    __safe_for_unpickling__ = True
    pass
