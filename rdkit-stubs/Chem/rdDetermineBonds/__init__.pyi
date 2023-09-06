"""Module containing a C++ implementation of the xyz2mol algorithm. This is based on xyz2mol: https://github.com/jensengroup/xyz2mol"""
from __future__ import annotations
import rdkit.Chem.rdDetermineBonds
import typing

__all__ = [
    "DetermineBondOrders",
    "DetermineBonds",
    "DetermineConnectivity"
]


def DetermineBondOrders( mol: Mol, charge: int = 0, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False) -> None:
    """
    DetermineBondOrders( mol: Mol, charge: int = 0, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False) -> None
        Assigns atomic connectivity to a molecule using atomic coordinates, 
        disregarding pre-existing bonds
        
        Args:
           mol : the molecule of interest; it must have a 3D conformer
           charge : (optional) the charge of the molecule; it must be provided if
               the Hueckel method is used and charge is non-zero
           allowChargedFragments : (optional) if this is true, formal charges
               will be placed on atoms according to their valency; otherwise, radical
               electrons will be placed on the atoms
           embedChiral : (optional) if this is true,
              chirality information will be embedded into the molecule; the function calls
              sanitizeMol() when this is true
           useAtomMap : (optional) if this is true, an atom map will be created for the 
              molecule
        

        C++ signature :
            void DetermineBondOrders(RDKit::ROMol {lvalue} [,int=0 [,bool=True [,bool=True [,bool=False]]]])
    """
def DetermineBonds( mol: Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False) -> None:
    """
    DetermineBonds( mol: Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False) -> None
        Assigns atomic connectivity to a molecule using atomic coordinates, 
        disregarding pre-existing bonds
        
        Args:
           mol : the molecule of interest; it must have a 3D conformer
           useHueckel : (optional) if this is  \c true, extended Hueckel theory
               will be used to determine connectivity rather than the van der Waals method
           charge : (optional) the charge of the molecule; it must be provided if
               the Hueckel method is used and charge is non-zero
           covFactor : (optional) the factor with which to multiply each covalent
               radius if the van der Waals method is used
           allowChargedFragments : (optional) if this is true, formal charges
               will be placed on atoms according to their valency; otherwise, radical
               electrons will be placed on the atoms
           embedChiral : (optional) if this is true,
              chirality information will be embedded into the molecule; the function calls
              sanitizeMol() when this is true
           useAtomMap : (optional) if this is true, an atom map will be created for the 
              molecule
        

        C++ signature :
            void DetermineBonds(RDKit::ROMol {lvalue} [,bool=False [,int=0 [,double=1.3 [,bool=True [,bool=True [,bool=False]]]]]])
    """
def DetermineConnectivity( mol: Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3) -> None:
    """
    DetermineConnectivity( mol: Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3) -> None
        Assigns atomic connectivity to a molecule using atomic coordinates, 
        disregarding pre-existing bonds
        
        Args:
           mol : the molecule of interest; it must have a 3D conformer
           useHueckel : (optional) if this is  \c true, extended Hueckel theory
               will be used to determine connectivity rather than the van der Waals method
           charge : (optional) the charge of the molecule; it must be provided if
               the Hueckel method is used and charge is non-zero
           covFactor : (optional) the factor with which to multiply each covalent
               radius if the van der Waals method is used
        

        C++ signature :
            void DetermineConnectivity(RDKit::ROMol {lvalue} [,bool=False [,int=0 [,double=1.3]]])
    """
