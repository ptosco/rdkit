"""Module containing classes and functions for enumerating molecules"""
from __future__ import annotations
import rdkit.Chem.rdMolEnumerator
import typing
import Boost.Python

__all__ = [
    "Enumerate",
    "EnumeratorType",
    "MolEnumeratorParams"
]


class EnumeratorType(Boost.Python.enum, int):
    LinkNode = rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode
    PositionVariation = rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation
    RepeatUnit = rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit
    __slots__ = ()
    names = {'LinkNode': rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode, 'PositionVariation': rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation, 'RepeatUnit': rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit}
    values = {0: rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode, 1: rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation, 2: rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit}
    pass
class MolEnumeratorParams(Boost.Python.instance):
    """
    Molecular enumerator parameters
    """
    @staticmethod
    def SetEnumerationOperator( arg1: MolEnumeratorParams, arg2: EnumeratorType) -> None: 
        """
        SetEnumerationOperator( arg1: MolEnumeratorParams, arg2: EnumeratorType) -> None
            set the operator to be used for enumeration

            C++ signature :
                void SetEnumerationOperator(RDKit::MolEnumerator::MolEnumeratorParams*,(anonymous namespace)::EnumeratorTypes)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void* __init__(boost::python::api::object,(anonymous namespace)::EnumeratorTypes)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, arg2: EnumeratorType) -> object: ...
    @property
    def doRandom(self) -> None:
        """
        do random enumeration (not yet implemented

        :type: None
        """
    @property
    def maxToEnumerate(self) -> None:
        """
        maximum number of molecules to enumerate

        :type: None
        """
    @property
    def randomSeed(self) -> None:
        """
        seed for the random enumeration (not yet implemented

        :type: None
        """
    @property
    def sanitize(self) -> None:
        """
        sanitize molecules after enumeration

        :type: None
        """
    __instance_size__ = 64
    pass
@typing.overload
def Enumerate( mol: Mol, maxPerOperation: int = 0) -> MolBundle:
    """
    Enumerate( mol: Mol, maxPerOperation: int = 0) -> MolBundle
        do an enumeration and return a MolBundle.
          If maxPerOperation is >0 that will be used as the maximum number of molecules which 
            can be returned by any given operation.
        Limitations:
          - the current implementation does not support molecules which include both
            SRUs and LINKNODEs
          - Overlapping SRUs, i.e. where one monomer is contained within another, are
            not supported

        C++ signature :
            RDKit::MolBundle* Enumerate(RDKit::ROMol [,unsigned int=0])

        C++ signature :
            RDKit::MolBundle* Enumerate(RDKit::ROMol,RDKit::MolEnumerator::MolEnumeratorParams)
    """
@typing.overload
def Enumerate( mol: Mol, enumParams: MolEnumeratorParams) -> MolBundle:
    pass
