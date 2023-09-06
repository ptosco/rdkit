"""Module containing rdFreeSASA classes and functions."""
from __future__ import annotations
import rdkit.Chem.rdFreeSASA
import typing
import Boost.Python

__all__ = [
    "APolar",
    "CalcSASA",
    "LeeRichards",
    "MakeFreeSasaAPolarAtomQuery",
    "MakeFreeSasaPolarAtomQuery",
    "NACCESS",
    "OONS",
    "Polar",
    "Protor",
    "SASAAlgorithm",
    "SASAClass",
    "SASAClassifier",
    "SASAOpts",
    "ShrakeRupley",
    "Unclassified",
    "classifyAtoms"
]


class SASAAlgorithm(Boost.Python.enum, int):
    LeeRichards = rdkit.Chem.rdFreeSASA.SASAAlgorithm.LeeRichards
    ShrakeRupley = rdkit.Chem.rdFreeSASA.SASAAlgorithm.ShrakeRupley
    __slots__ = ()
    names = {'LeeRichards': rdkit.Chem.rdFreeSASA.SASAAlgorithm.LeeRichards, 'ShrakeRupley': rdkit.Chem.rdFreeSASA.SASAAlgorithm.ShrakeRupley}
    values = {0: rdkit.Chem.rdFreeSASA.SASAAlgorithm.LeeRichards, 1: rdkit.Chem.rdFreeSASA.SASAAlgorithm.ShrakeRupley}
    pass
class SASAClass(Boost.Python.enum, int):
    APolar = rdkit.Chem.rdFreeSASA.SASAClass.APolar
    Polar = rdkit.Chem.rdFreeSASA.SASAClass.Polar
    Unclassified = rdkit.Chem.rdFreeSASA.SASAClass.Unclassified
    __slots__ = ()
    names = {'Unclassified': rdkit.Chem.rdFreeSASA.SASAClass.Unclassified, 'APolar': rdkit.Chem.rdFreeSASA.SASAClass.APolar, 'Polar': rdkit.Chem.rdFreeSASA.SASAClass.Polar}
    values = {0: rdkit.Chem.rdFreeSASA.SASAClass.Unclassified, 1: rdkit.Chem.rdFreeSASA.SASAClass.APolar, 2: rdkit.Chem.rdFreeSASA.SASAClass.Polar}
    pass
class SASAClassifier(Boost.Python.enum, int):
    NACCESS = rdkit.Chem.rdFreeSASA.SASAClassifier.NACCESS
    OONS = rdkit.Chem.rdFreeSASA.SASAClassifier.OONS
    Protor = rdkit.Chem.rdFreeSASA.SASAClassifier.Protor
    __slots__ = ()
    names = {'Protor': rdkit.Chem.rdFreeSASA.SASAClassifier.Protor, 'NACCESS': rdkit.Chem.rdFreeSASA.SASAClassifier.NACCESS, 'OONS': rdkit.Chem.rdFreeSASA.SASAClassifier.OONS}
    values = {0: rdkit.Chem.rdFreeSASA.SASAClassifier.Protor, 1: rdkit.Chem.rdFreeSASA.SASAClassifier.NACCESS, 2: rdkit.Chem.rdFreeSASA.SASAClassifier.OONS}
    pass
class SASAOpts(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None
            Constructor takes no arguments

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,FreeSASA::SASAOpts::Algorithm,FreeSASA::SASAOpts::Classifier)

            C++ signature :
                void __init__(_object*,FreeSASA::SASAOpts::Algorithm,FreeSASA::SASAOpts::Classifier,double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: SASAAlgorithm, arg3: SASAClassifier) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: SASAAlgorithm, arg3: SASAClassifier, arg4: float) -> None: ...
    @property
    def algorithm(self) -> None:
        """
        :type: None
        """
    @property
    def classifier(self) -> None:
        """
        :type: None
        """
    @property
    def probeRadius(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 40
    pass
def CalcSASA( mol: Mol, radii: AtomPairsParameters, confIdx: int = -1, query: Atom = None, opts: SASAOpts = SASAOpts()) -> float:
    """
    CalcSASA( mol: Mol, radii: AtomPairsParameters, confIdx: int = -1, query: Atom = None, opts: SASAOpts = SASAOpts()) -> float
        Compute the Solvent Accessible Surface Area using the FreeSASA library
        ARGUMENTS:
          - mol: The molecule to compute.
          - radii:  A list of atom raddii where radii[atom.GetIdx()] is the radius of the atom
                    These can be passed in or calculated with classifyAtoms for some proteins
          - confIdx: Specify the conformer to use for the 3D geometry  [default -1]
          - query: Pass along a query atom to compute the SASA for a subset of atoms.
                   precanned query atoms can be made with MakeFreeSasaPolarAtomQuery and
                   MakeFreeSasaAPolarAtomQuery for classified polar and apolar atoms respectively.
          - opts: a SASAOpts class specifying the algorithm to use
        
        RETURNS:
        The computed solvent accessible surface area.
        

        C++ signature :
            double CalcSASA(RDKit::ROMol,boost::python::api::object [,int=-1 [,RDKit::Atom const*=None [,FreeSASA::SASAOpts=<rdkit.Chem.rdFreeSASA.SASAOpts object at 0x7fd051ac4f90>]]])
    """
def MakeFreeSasaAPolarAtomQuery() -> QueryAtom:
    """
    MakeFreeSasaAPolarAtomQuery() -> QueryAtom
        Returns an APolar atom query for use with CalcSASA.  An apolar atom has the SASAClass
        and SASAClassName set to the APOLAR class.  (see classifyAtoms)

        C++ signature :
            RDKit::QueryAtom const* MakeFreeSasaAPolarAtomQuery()
    """
def MakeFreeSasaPolarAtomQuery() -> QueryAtom:
    """
    MakeFreeSasaPolarAtomQuery() -> QueryAtom
        Returns a polar atom query for use with CalcSASA.  An polar atom has the SASAClass
        and SASAClassName set to the POLAR class.  (see classifyAtoms)

        C++ signature :
            RDKit::QueryAtom const* MakeFreeSasaPolarAtomQuery()
    """
def classifyAtoms( mol: Mol, options: SASAOpts = SASAOpts()) -> object:
    """
    classifyAtoms( mol: Mol, options: SASAOpts = SASAOpts()) -> object
        Classify the atoms in the molecule returning their radii if possible.
        ARGUMENTS:
           - mol: molecule to classify
           - options: FreeSASA options class specifying the classification method.
                       Current classifiers are Protor, NACCESS and OONS
                       classification is stored as atom property 'SASAClass' for the integer value
                        and 'SASAClassName' for the string name of the class, Polar, APolar...
        
        RETURNS:
          list of radii where radii[atom.GetIdx()] is the radii of the atom.
          If classification fails, NONE is returned
        

        C++ signature :
            boost::python::api::object classifyAtoms(RDKit::ROMol {lvalue} [,FreeSASA::SASAOpts=<rdkit.Chem.rdFreeSASA.SASAOpts object at 0x7fd051ac4dd0>])
    """
APolar = rdkit.Chem.rdFreeSASA.SASAClass.APolar
LeeRichards = rdkit.Chem.rdFreeSASA.SASAAlgorithm.LeeRichards
NACCESS = rdkit.Chem.rdFreeSASA.SASAClassifier.NACCESS
OONS = rdkit.Chem.rdFreeSASA.SASAClassifier.OONS
Polar = rdkit.Chem.rdFreeSASA.SASAClass.Polar
Protor = rdkit.Chem.rdFreeSASA.SASAClassifier.Protor
ShrakeRupley = rdkit.Chem.rdFreeSASA.SASAAlgorithm.ShrakeRupley
Unclassified = rdkit.Chem.rdFreeSASA.SASAClass.Unclassified
