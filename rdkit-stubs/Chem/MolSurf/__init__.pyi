""" Exposes functionality for MOE-like approximate molecular surface area
descriptors.

  The MOE-like VSA descriptors are also calculated here

"""
from __future__ import annotations
import rdkit.Chem.MolSurf
import typing
import bisect
import numpy
import rdkit.Chem
import rdkit.Chem.Crippen
import rdkit.Chem.rdMolDescriptors
import rdkit.Chem.rdPartialCharges
import rdkit.Chem.rdchem
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Chem",
    "Crippen",
    "LabuteASA",
    "PEOE_VSA1",
    "PEOE_VSA10",
    "PEOE_VSA11",
    "PEOE_VSA12",
    "PEOE_VSA13",
    "PEOE_VSA14",
    "PEOE_VSA2",
    "PEOE_VSA3",
    "PEOE_VSA4",
    "PEOE_VSA5",
    "PEOE_VSA6",
    "PEOE_VSA7",
    "PEOE_VSA8",
    "PEOE_VSA9",
    "PEOE_VSA_",
    "SMR_VSA1",
    "SMR_VSA10",
    "SMR_VSA2",
    "SMR_VSA3",
    "SMR_VSA4",
    "SMR_VSA5",
    "SMR_VSA6",
    "SMR_VSA7",
    "SMR_VSA8",
    "SMR_VSA9",
    "SMR_VSA_",
    "SlogP_VSA1",
    "SlogP_VSA10",
    "SlogP_VSA11",
    "SlogP_VSA12",
    "SlogP_VSA2",
    "SlogP_VSA3",
    "SlogP_VSA4",
    "SlogP_VSA5",
    "SlogP_VSA6",
    "SlogP_VSA7",
    "SlogP_VSA8",
    "SlogP_VSA9",
    "SlogP_VSA_",
    "TPSA",
    "bisect",
    "bondScaleFacts",
    "chgBins",
    "logpBins",
    "mrBins",
    "numpy",
    "ptable",
    "pyLabuteASA",
    "pyPEOE_VSA_",
    "pySMR_VSA_",
    "pySlogP_VSA_",
    "rdMolDescriptors",
    "rdPartialCharges"
]


def PEOE_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    PEOE_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list PEOE_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
def SMR_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    SMR_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list SMR_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
def SlogP_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    SlogP_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list SlogP_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
bondScaleFacts = [0.1, 0, 0.2, 0.3]
chgBins = [-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
logpBins = [-0.4, -0.2, 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
mrBins = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0]
ptable: rdkit.Chem.rdchem.PeriodicTable
LabuteASA = rdkit.Chem.MolSurf.<lambda>
PEOE_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA11 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA12 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA13 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA14 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA11 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA12 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
TPSA = rdkit.Chem.MolSurf.<lambda>
