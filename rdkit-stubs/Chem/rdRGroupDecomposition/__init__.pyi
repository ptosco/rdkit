"""Module containing RGroupDecomposition classes and functions."""
from __future__ import annotations
import rdkit.Chem.rdRGroupDecomposition
import typing
import Boost.Python

__all__ = [
    "AtomIndexLabels",
    "AtomMap",
    "AtomMapLabels",
    "AutoDetect",
    "DummyAtomLabels",
    "Exhaustive",
    "FingerprintVariance",
    "GA",
    "Greedy",
    "GreedyChunks",
    "Isotope",
    "IsotopeLabels",
    "MCS",
    "MDLRGroup",
    "MDLRGroupLabels",
    "Match",
    "NoAlignment",
    "NoSymmetrization",
    "None",
    "RGroupCoreAlignment",
    "RGroupDecompose",
    "RGroupDecomposition",
    "RGroupDecompositionParameters",
    "RGroupLabelling",
    "RGroupLabels",
    "RGroupMatching",
    "RGroupScore",
    "RelabelDuplicateLabels"
]


class RGroupCoreAlignment(Boost.Python.enum, int):
    MCS = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS
    NoAlignment = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment
    None = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None
    __slots__ = ()
    names = {'None': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None, 'NoAlignment': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment, 'MCS': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS}
    values = {0: rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment, 1: rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS}
    pass
class RGroupDecomposition(Boost.Python.instance):
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
      OPTIONS:
        - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCS
                               If set to MCS, cores labels are mapped to each other using their
                               Maximum common substructure overlap.
        - RGroupLabels: optionally set where the rgroup labels to use are encoded.
                         RGroupLabels.IsotopeLabels - labels are stored on isotopes
                         RGroupLabels.AtomMapLabels - labels are stored on atommaps
                         RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
                         RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
                         RGroupLabels.AtomIndexLabels - use the atom index as the label
                         RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
                         RGroupLabels.AutoDetect - auto detect the label [default]
           Note: in all cases, any rgroups found on unlabelled atoms will be automatically
                  labelled.
        - RGroupLabelling: choose where the rlabels are stored on the decomposition
                            RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
                            RGroupLabelling.Isotope - store rgroups on the isotope
                            RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)
                           default: AtomMap | MDLRGroup
        - onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
        - removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
        - removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
        - removeHydrogensPostMatch: remove all hydrogens from the output molecules
        - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more
        - doTautomers: match all tautomers of a core against each input structure
        - doEnumeration: expand input cores into enumerated mol bundles
        -allowMultipleRGroupsOnUnlabelled: permit more that one rgroup to be attached to an unlabelled core atom
    """
    @staticmethod
    def Add( arg1: RGroupDecomposition, arg2: Mol) -> int: 
        """
        Add( arg1: RGroupDecomposition, arg2: Mol) -> int

            C++ signature :
                int Add(RDKit::RGroupDecompositionHelper {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def GetRGroupLabels( arg1: RGroupDecomposition) -> list: 
        """
        GetRGroupLabels( arg1: RGroupDecomposition) -> list
            Return the current list of found rgroups.
            Note, Process() should be called first

            C++ signature :
                boost::python::list GetRGroupLabels(RDKit::RGroupDecompositionHelper {lvalue})
        """
    @staticmethod
    def GetRGroupsAsColumns( arg1: RGroupDecomposition, asSmiles: bool = False) -> dict: 
        """
        GetRGroupsAsColumns( arg1: RGroupDecomposition, asSmiles: bool = False) -> dict
            Return the rgroups as columns (note: can be fed directrly into a pandas datatable)
              ARGUMENTS:
               - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
                Column structure:
                   columns[rgroup_label] = [ mols_or_smiles ]
            

            C++ signature :
                boost::python::dict GetRGroupsAsColumns(RDKit::RGroupDecompositionHelper {lvalue} [,bool=False])
        """
    @staticmethod
    def GetRGroupsAsRows( arg1: RGroupDecomposition, asSmiles: bool = False) -> list: 
        """
        GetRGroupsAsRows( arg1: RGroupDecomposition, asSmiles: bool = False) -> list
            Return the rgroups as rows (note: can be fed directrly into a pandas datatable)
              ARGUMENTS:
               - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
                Row structure:
                   rows[idx] = {rgroup_label: molecule_or_smiles}
            

            C++ signature :
                boost::python::list GetRGroupsAsRows(RDKit::RGroupDecompositionHelper {lvalue} [,bool=False])
        """
    @staticmethod
    def Process( arg1: RGroupDecomposition) -> bool: 
        """
        Process( arg1: RGroupDecomposition) -> bool
            Process the rgroups (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)

            C++ signature :
                bool Process(RDKit::RGroupDecompositionHelper {lvalue})
        """
    @staticmethod
    def ProcessAndScore( arg1: RGroupDecomposition) -> tuple: 
        """
        ProcessAndScore( arg1: RGroupDecomposition) -> tuple
            Process the rgroups and returns the score (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)

            C++ signature :
                boost::python::tuple ProcessAndScore(RDKit::RGroupDecompositionHelper {lvalue})
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: AtomPairsParameters) -> None: 
        """
        __init__( arg1: object, arg2: AtomPairsParameters) -> None
            Construct from a molecule or sequence of molecules

            C++ signature :
                void __init__(_object*,boost::python::api::object)

            C++ signature :
                void __init__(_object*,boost::python::api::object,RDKit::RGroupDecompositionParameters)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: AtomPairsParameters, arg3: RGroupDecompositionParameters) -> None: ...
    __instance_size__ = 32
    pass
class RGroupDecompositionParameters(Boost.Python.instance):
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
      OPTIONS:
        - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCS
                               If set to MCS, cores labels are mapped to each other using their
                               Maximum common substructure overlap.
        - RGroupLabels: optionally set where the rgroup labels to use are encoded.
                         RGroupLabels.IsotopeLabels - labels are stored on isotopes
                         RGroupLabels.AtomMapLabels - labels are stored on atommaps
                         RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
                         RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
                         RGroupLabels.AtomIndexLabels - use the atom index as the label
                         RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
                         RGroupLabels.AutoDetect - auto detect the label [default]
           Note: in all cases, any rgroups found on unlabelled atoms will be automatically
                  labelled.
        - RGroupLabelling: choose where the rlabels are stored on the decomposition
                            RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
                            RGroupLabelling.Isotope - store rgroups on the isotope
                            RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)
                           default: AtomMap | MDLRGroup
        - onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
        - removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
        - removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
        - removeHydrogensPostMatch: remove all hydrogens from the output molecules
        - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more
        - doTautomers: match all tautomers of a core against each input structure
        - doEnumeration: expand input cores into enumerated mol bundles
        -allowMultipleRGroupsOnUnlabelled: permit more that one rgroup to be attached to an unlabelled core atom
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None
            Constructor, takes no arguments

            C++ signature :
                void __init__(_object*)
        """
    @property
    def alignment(self) -> None:
        """
        :type: None
        """
    @property
    def allowMultipleRGroupsOnUnlabelled(self) -> None:
        """
        :type: None
        """
    @property
    def allowNonTerminalRGroups(self) -> None:
        """
        :type: None
        """
    @property
    def chunkSize(self) -> None:
        """
        :type: None
        """
    @property
    def doEnumeration(self) -> None:
        """
        :type: None
        """
    @property
    def doTautomers(self) -> None:
        """
        :type: None
        """
    @property
    def gaMaximumOperations(self) -> None:
        """
        :type: None
        """
    @property
    def gaNumberOperationsWithoutImprovement(self) -> None:
        """
        :type: None
        """
    @property
    def gaNumberRuns(self) -> None:
        """
        :type: None
        """
    @property
    def gaParallelRuns(self) -> None:
        """
        :type: None
        """
    @property
    def gaPopulationSize(self) -> None:
        """
        :type: None
        """
    @property
    def gaRandomSeed(self) -> None:
        """
        :type: None
        """
    @property
    def labels(self) -> None:
        """
        :type: None
        """
    @property
    def matchingStrategy(self) -> None:
        """
        :type: None
        """
    @property
    def onlyMatchAtRGroups(self) -> None:
        """
        :type: None
        """
    @property
    def removeAllHydrogenRGroups(self) -> None:
        """
        :type: None
        """
    @property
    def removeAllHydrogenRGroupsAndLabels(self) -> None:
        """
        :type: None
        """
    @property
    def removeHydrogensPostMatch(self) -> None:
        """
        :type: None
        """
    @property
    def rgroupLabelling(self) -> None:
        """
        :type: None
        """
    @property
    def scoreMethod(self) -> None:
        """
        :type: None
        """
    @property
    def substructMatchParams(self) -> None:
        """
        :type: None
        """
    @property
    def timeout(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 192
    pass
class RGroupLabelling(Boost.Python.enum, int):
    AtomMap = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap
    Isotope = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope
    MDLRGroup = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup
    __slots__ = ()
    names = {'AtomMap': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap, 'Isotope': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope, 'MDLRGroup': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup}
    values = {1: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap, 2: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope, 4: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup}
    pass
class RGroupLabels(Boost.Python.enum, int):
    AtomIndexLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels
    AtomMapLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels
    AutoDetect = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect
    DummyAtomLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels
    IsotopeLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels
    MDLRGroupLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels
    RelabelDuplicateLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels
    __slots__ = ()
    names = {'IsotopeLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels, 'AtomMapLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels, 'AtomIndexLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels, 'RelabelDuplicateLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels, 'MDLRGroupLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels, 'DummyAtomLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels, 'AutoDetect': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect}
    values = {1: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels, 2: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels, 4: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels, 8: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels, 16: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels, 32: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels, 255: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect}
    pass
class RGroupMatching(Boost.Python.enum, int):
    Exhaustive = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive
    GA = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA
    Greedy = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy
    GreedyChunks = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks
    NoSymmetrization = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization
    __slots__ = ()
    names = {'Greedy': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy, 'GreedyChunks': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks, 'Exhaustive': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive, 'NoSymmetrization': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization, 'GA': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA}
    values = {1: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy, 2: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks, 4: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive, 8: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization, 16: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA}
    pass
class RGroupScore(Boost.Python.enum, int):
    FingerprintVariance = rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance
    Match = rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match
    __slots__ = ()
    names = {'Match': rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match, 'FingerprintVariance': rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance}
    values = {1: rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match, 4: rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance}
    pass
def RGroupDecompose( cores: AtomPairsParameters, mols: AtomPairsParameters, asSmiles: bool = False, asRows: bool = True, options: RGroupDecompositionParameters = RGroupDecompositionParameters()) -> object:
    """
    RGroupDecompose( cores: AtomPairsParameters, mols: AtomPairsParameters, asSmiles: bool = False, asRows: bool = True, options: RGroupDecompositionParameters = RGroupDecompositionParameters()) -> object
        Decompose a collecion of molecules into their Rgroups
          ARGUMENTS:
            - cores: a set of cores from most to least specific.
                     See RGroupDecompositionParameters for more details
                     on how the cores can be labelled
            - mols: the molecules to be decomposed
            - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
            - asRows: return the results as rows (default) otherwise return columns
        
          RETURNS: row_or_column_results, unmatched
        
            Row structure:
               rows[idx] = {rgroup_label: molecule_or_smiles}
            Column structure:
               columns[rgroup_label] = [ mols_or_smiles ]
        
            unmatched is a vector of indices in the input mols that were not matched.
        

        C++ signature :
            boost::python::api::object RGroupDecompose(boost::python::api::object,boost::python::api::object [,bool=False [,bool=True [,RDKit::RGroupDecompositionParameters=<rdkit.Chem.rdRGroupDecomposition.RGroupDecompositionParameters object at 0x7f2c611be040>]]])
    """
AtomIndexLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels
AtomMap = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap
AtomMapLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels
AutoDetect = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect
DummyAtomLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels
Exhaustive = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive
FingerprintVariance = rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance
GA = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA
Greedy = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy
GreedyChunks = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks
Isotope = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope
IsotopeLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels
MCS = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS
MDLRGroup = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup
MDLRGroupLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels
Match = rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match
NoAlignment = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment
NoSymmetrization = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization
None = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None
RelabelDuplicateLabels = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels
