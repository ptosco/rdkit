"""Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next."""
from __future__ import annotations
import rdkit.Chem.rdEHTTools
import typing
import Boost.Python

__all__ = [
    "EHTResults",
    "RunMol"
]


class EHTResults(Boost.Python.instance):
    @staticmethod
    def GetAtomicCharges( arg1: EHTResults) -> object: 
        """
        GetAtomicCharges( arg1: EHTResults) -> object
            returns the calculated atomic charges

            C++ signature :
                _object* GetAtomicCharges(RDKit::EHTTools::EHTResults {lvalue})
        """
    @staticmethod
    def GetHamiltonian( arg1: EHTResults) -> object: 
        """
        GetHamiltonian( arg1: EHTResults) -> object
            returns the symmetric Hamiltonian matrix

            C++ signature :
                _object* GetHamiltonian(RDKit::EHTTools::EHTResults {lvalue})
        """
    @staticmethod
    def GetOrbitalEnergies( arg1: EHTResults) -> object: 
        """
        GetOrbitalEnergies( arg1: EHTResults) -> object
            returns the energies of the molecular orbitals as a vector

            C++ signature :
                _object* GetOrbitalEnergies(RDKit::EHTTools::EHTResults {lvalue})
        """
    @staticmethod
    def GetOverlapMatrix( arg1: EHTResults) -> object: 
        """
        GetOverlapMatrix( arg1: EHTResults) -> object
            returns the symmetric overlap matrix

            C++ signature :
                _object* GetOverlapMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
    @staticmethod
    def GetReducedChargeMatrix( arg1: EHTResults) -> object: 
        """
        GetReducedChargeMatrix( arg1: EHTResults) -> object
            returns the reduced charge matrix

            C++ signature :
                _object* GetReducedChargeMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
    @staticmethod
    def GetReducedOverlapPopulationMatrix( arg1: EHTResults) -> object: 
        """
        GetReducedOverlapPopulationMatrix( arg1: EHTResults) -> object
            returns the reduced overlap population matrix

            C++ signature :
                _object* GetReducedOverlapPopulationMatrix(RDKit::EHTTools::EHTResults {lvalue})
        """
    @property
    def fermiEnergy(self) -> None:
        """
        :type: None
        """
    @property
    def numElectrons(self) -> None:
        """
        :type: None
        """
    @property
    def numOrbitals(self) -> None:
        """
        :type: None
        """
    @property
    def totalEnergy(self) -> None:
        """
        :type: None
        """
    pass
def RunMol( mol: Mol, confId: int = -1, keepOverlapAndHamiltonianMatrices: bool = False) -> tuple:
    """
    RunMol( mol: Mol, confId: int = -1, keepOverlapAndHamiltonianMatrices: bool = False) -> tuple
        Runs an extended Hueckel calculation for a molecule.
        The molecule should have at least one conformation
        
        ARGUMENTS:
           - mol: molecule to use
           - confId: (optional) conformation to use
           - keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap 
             and hamiltonian matrices in the EHTResults object
        
        RETURNS: a 2-tuple:
           - a boolean indicating whether or not the calculation succeeded
           - an EHTResults object with the results
        

        C++ signature :
            boost::python::tuple RunMol(RDKit::ROMol [,int=-1 [,bool=False]])
    """
