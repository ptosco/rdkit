""" A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation

"""
from __future__ import annotations
import rdkit.Chem
import typing
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.Chem.rdmolops import AdjustQueryWhichFlags
from rdkit.Chem.rdmolops import AromaticityModel
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.rdchem import AtomKekulizeException
from rdkit.Chem.rdchem import AtomMonomerInfo
from rdkit.Chem.rdchem import AtomMonomerType
from rdkit.Chem.rdchem import AtomPDBResidueInfo
from rdkit.Chem.rdchem import AtomSanitizeException
from rdkit.Chem.rdchem import AtomValenceException
from rdkit.Chem.rdchem import Bond
from rdkit.Chem.rdchem import BondDir
from rdkit.Chem.rdchem import BondStereo
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolops import BondWedgingParameters
from rdkit.Chem.rdmolfiles import CXSmilesFields
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdchem import CompositeQueryType
from rdkit.Chem.rdchem import Conformer
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdchem import FixedMolSizeMolBundle
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.inchi import InchiReadWriteError
from rdkit.Chem.rdMolInterchange import JSONParseParameters
from rdkit.Chem.rdMolInterchange import JSONWriteParameters
from rdkit.Chem.rdchem import KekulizeException
from rdkit.Chem.rdmolfiles import MaeMolSupplier
from rdkit.Chem.rdmolfiles import MaeWriter
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdchem import MolBundle
from rdkit.Chem.rdchem import MolSanitizeException
from rdkit.Chem.rdmolops import MolzipLabel
from rdkit.Chem.rdmolops import MolzipParams
from rdkit.Chem.rdmolfiles import MultithreadedSDMolSupplier
from rdkit.Chem.rdmolfiles import MultithreadedSmilesMolSupplier
from rdkit.Chem.rdmolfiles import PDBWriter
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem.rdchem import PropertyPickleOptions
from rdkit.Chem.rdchem import QueryAtom
from rdkit.Chem.rdchem import QueryBond
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdmolops import RemoveHsParameters
from rdkit.Chem.rdchem import ResonanceFlags
from rdkit.Chem.rdchem import ResonanceMolSupplier
from rdkit.Chem.rdchem import ResonanceMolSupplierCallback
from rdkit.Chem.rdchem import RingInfo
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem.rdmolops import SanitizeFlags
from rdkit.Chem.rdmolfiles import SmartsParserParams
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from rdkit.Chem.rdmolfiles import SmilesParserParams
from rdkit.Chem.rdmolfiles import SmilesWriteParams
from rdkit.Chem.rdmolfiles import SmilesWriter
from rdkit.Chem.rdmolops import StereoBondThresholds
from rdkit.Chem.rdchem import StereoDescriptor
from rdkit.Chem.rdchem import StereoGroup
from rdkit.Chem.rdchem import StereoGroupType
from rdkit.Chem.rdchem import StereoGroup_vect
from rdkit.Chem.rdchem import StereoInfo
from rdkit.Chem.rdchem import StereoSpecified
from rdkit.Chem.rdchem import StereoType
from rdkit.Chem.rdchem import SubstanceGroup
from rdkit.Chem.rdchem import SubstanceGroupAttach
from rdkit.Chem.rdchem import SubstanceGroupCState
from rdkit.Chem.rdchem import SubstanceGroup_VECT
from rdkit.Chem.rdchem import SubstructMatchParameters
from rdkit.Chem.rdmolfiles import TDTMolSupplier
from rdkit.Chem.rdmolfiles import TDTWriter
import rdkit.Chem.rdchem
import rdkit.Chem.rdmolops
import rdkit.DataStructs
import rdkit.Geometry.rdGeometry
import rdkit.RDConfig
import rdkit.rdBase

__all__ = [
    "ADJUST_IGNOREALL",
    "ADJUST_IGNORECHAINS",
    "ADJUST_IGNOREDUMMIES",
    "ADJUST_IGNOREMAPPED",
    "ADJUST_IGNORENONDUMMIES",
    "ADJUST_IGNORENONE",
    "ADJUST_IGNORERINGS",
    "ALLOW_CHARGE_SEPARATION",
    "ALLOW_INCOMPLETE_OCTETS",
    "AROMATICITY_CUSTOM",
    "AROMATICITY_DEFAULT",
    "AROMATICITY_MDL",
    "AROMATICITY_RDKIT",
    "AROMATICITY_SIMPLE",
    "AddHs",
    "AddMetadataToPNGFile",
    "AddMetadataToPNGString",
    "AddMolSubstanceGroup",
    "AddRecursiveQuery",
    "AddWavyBondsForStereoAny",
    "AdjustQueryParameters",
    "AdjustQueryProperties",
    "AdjustQueryPropertiesWithGenericGroups",
    "AdjustQueryWhichFlags",
    "AllProps",
    "AromaticityModel",
    "AssignAtomChiralTagsFromMolParity",
    "AssignAtomChiralTagsFromStructure",
    "AssignCIPLabels",
    "AssignChiralTypesFromBondDirs",
    "AssignRadicals",
    "AssignStereochemistry",
    "AssignStereochemistryFrom3D",
    "Atom",
    "AtomFromSmarts",
    "AtomFromSmiles",
    "AtomKekulizeException",
    "AtomMonomerInfo",
    "AtomMonomerType",
    "AtomPDBResidueInfo",
    "AtomProps",
    "AtomSanitizeException",
    "AtomValenceException",
    "Bond",
    "BondDir",
    "BondFromSmarts",
    "BondFromSmiles",
    "BondProps",
    "BondStereo",
    "BondType",
    "BondWedgingParameters",
    "CHI_ALLENE",
    "CHI_OCTAHEDRAL",
    "CHI_OTHER",
    "CHI_SQUAREPLANAR",
    "CHI_TETRAHEDRAL",
    "CHI_TETRAHEDRAL_CCW",
    "CHI_TETRAHEDRAL_CW",
    "CHI_TRIGONALBIPYRAMIDAL",
    "CHI_UNSPECIFIED",
    "COMPOSITE_AND",
    "COMPOSITE_OR",
    "COMPOSITE_XOR",
    "CXSmilesFields",
    "CanonSmiles",
    "CanonicalRankAtoms",
    "CanonicalRankAtomsInFragment",
    "CanonicalizeEnhancedStereo",
    "ChiralType",
    "Cleanup",
    "CleanupOrganometallics",
    "ClearMolSubstanceGroups",
    "CombineMols",
    "CompositeQueryType",
    "ComputedProps",
    "Conformer",
    "ConvertGenericQueriesToSubstanceGroups",
    "CoordsAsDouble",
    "CreateAtomBoolPropertyList",
    "CreateAtomDoublePropertyList",
    "CreateAtomIntPropertyList",
    "CreateAtomStringPropertyList",
    "CreateMolDataSubstanceGroup",
    "CreateMolSubstanceGroup",
    "CreateStereoGroup",
    "DataStructs",
    "DativeBondsToHaptic",
    "DeleteSubstructs",
    "DetectBondStereoChemistry",
    "DetectBondStereochemistry",
    "DetectChemistryProblems",
    "EditableMol",
    "FastFindRings",
    "FindAllPathsOfLengthN",
    "FindAllSubgraphsOfLengthMToN",
    "FindAllSubgraphsOfLengthN",
    "FindAtomEnvironmentOfRadiusN",
    "FindMolChiralCenters",
    "FindPotentialStereo",
    "FindPotentialStereoBonds",
    "FindRingFamilies",
    "FindUniqueSubgraphsOfLengthN",
    "FixedMolSizeMolBundle",
    "ForwardSDMolSupplier",
    "ForwardStereoGroupIds",
    "FragmentOnBRICSBonds",
    "FragmentOnBonds",
    "FragmentOnSomeBonds",
    "Get3DDistanceMatrix",
    "GetAdjacencyMatrix",
    "GetAllowNontetrahedralChirality",
    "GetAtomAlias",
    "GetAtomRLabel",
    "GetAtomValue",
    "GetDefaultPickleProperties",
    "GetDistanceMatrix",
    "GetFormalCharge",
    "GetMolFrags",
    "GetMolSubstanceGroupWithIdx",
    "GetMolSubstanceGroups",
    "GetMostSubstitutedCoreMatch",
    "GetPeriodicTable",
    "GetSSSR",
    "GetShortestPath",
    "GetSupplementalSmilesLabel",
    "GetSymmSSSR",
    "GetUseLegacyStereoPerception",
    "HapticBondsToDative",
    "HybridizationType",
    "INCHI_AVAILABLE",
    "InchiReadWriteError",
    "InchiToInchiKey",
    "JSONParseParameters",
    "JSONToMols",
    "JSONWriteParameters",
    "KEKULE_ALL",
    "Kekulize",
    "KekulizeException",
    "KekulizeIfPossible",
    "LayeredFingerprint",
    "LayeredFingerprint_substructLayers",
    "MaeMolSupplier",
    "MaeWriter",
    "MergeQueryHs",
    "MetadataFromPNGFile",
    "MetadataFromPNGString",
    "Mol",
    "MolAddRecursiveQueries",
    "MolBlockToInchi",
    "MolBlockToInchiAndAuxInfo",
    "MolBundle",
    "MolBundleCanSerialize",
    "MolFragmentToCXSmarts",
    "MolFragmentToCXSmiles",
    "MolFragmentToSmarts",
    "MolFragmentToSmiles",
    "MolFromFASTA",
    "MolFromHELM",
    "MolFromInchi",
    "MolFromMol2Block",
    "MolFromMol2File",
    "MolFromMolBlock",
    "MolFromMolFile",
    "MolFromMrvBlock",
    "MolFromMrvFile",
    "MolFromPDBBlock",
    "MolFromPDBFile",
    "MolFromPNGFile",
    "MolFromPNGString",
    "MolFromRDKitSVG",
    "MolFromSequence",
    "MolFromSmarts",
    "MolFromSmiles",
    "MolFromTPLBlock",
    "MolFromTPLFile",
    "MolFromXYZBlock",
    "MolFromXYZFile",
    "MolMetadataToPNGFile",
    "MolMetadataToPNGString",
    "MolProps",
    "MolSanitizeException",
    "MolToCMLBlock",
    "MolToCMLFile",
    "MolToCXSmarts",
    "MolToCXSmiles",
    "MolToFASTA",
    "MolToHELM",
    "MolToInchi",
    "MolToInchiAndAuxInfo",
    "MolToInchiKey",
    "MolToJSON",
    "MolToMolBlock",
    "MolToMolFile",
    "MolToMrvBlock",
    "MolToMrvFile",
    "MolToPDBBlock",
    "MolToPDBFile",
    "MolToRandomSmilesVect",
    "MolToSequence",
    "MolToSmarts",
    "MolToSmiles",
    "MolToTPLBlock",
    "MolToTPLFile",
    "MolToV3KMolBlock",
    "MolToV3KMolFile",
    "MolToXYZBlock",
    "MolToXYZFile",
    "MolsFromCDXML",
    "MolsFromCDXMLFile",
    "MolsFromPNGFile",
    "MolsFromPNGString",
    "MolsToJSON",
    "MolzipLabel",
    "MolzipParams",
    "MultithreadedSDMolSupplier",
    "MultithreadedSmilesMolSupplier",
    "MurckoDecompose",
    "NoConformers",
    "NoProps",
    "PDBWriter",
    "ParseMolQueryDefFile",
    "PathToSubmol",
    "PatternFingerprint",
    "PeriodicTable",
    "PrivateProps",
    "PropertyPickleOptions",
    "QueryAtom",
    "QueryAtomData",
    "QueryBond",
    "QuickSmartsMatch",
    "RDConfig",
    "RDKFingerprint",
    "RWMol",
    "ReapplyMolBlockWedging",
    "RemoveAllHs",
    "RemoveHs",
    "RemoveHsParameters",
    "RemoveStereochemistry",
    "RenumberAtoms",
    "ReplaceCore",
    "ReplaceSidechains",
    "ReplaceSubstructs",
    "ResonanceFlags",
    "ResonanceMolSupplier",
    "ResonanceMolSupplierCallback",
    "RingInfo",
    "SANITIZE_ADJUSTHS",
    "SANITIZE_ALL",
    "SANITIZE_CLEANUP",
    "SANITIZE_CLEANUPCHIRALITY",
    "SANITIZE_CLEANUP_ORGANOMETALLICS",
    "SANITIZE_FINDRADICALS",
    "SANITIZE_KEKULIZE",
    "SANITIZE_NONE",
    "SANITIZE_PROPERTIES",
    "SANITIZE_SETAROMATICITY",
    "SANITIZE_SETCONJUGATION",
    "SANITIZE_SETHYBRIDIZATION",
    "SANITIZE_SYMMRINGS",
    "SDMolSupplier",
    "SDWriter",
    "STEREO_ABSOLUTE",
    "STEREO_AND",
    "STEREO_OR",
    "SanitizeFlags",
    "SanitizeMol",
    "SetAllowNontetrahedralChirality",
    "SetAromaticity",
    "SetAtomAlias",
    "SetAtomRLabel",
    "SetAtomValue",
    "SetBondStereoFromDirections",
    "SetConjugation",
    "SetDefaultPickleProperties",
    "SetDoubleBondNeighborDirections",
    "SetGenericQueriesFromProperties",
    "SetHybridization",
    "SetSupplementalSmilesLabel",
    "SetTerminalAtomCoords",
    "SetUseLegacyStereoPerception",
    "SmartsParserParams",
    "SmilesMolSupplier",
    "SmilesMolSupplierFromText",
    "SmilesParserParams",
    "SmilesWriteParams",
    "SmilesWriter",
    "SortMatchesByDegreeOfCoreSubstitution",
    "SplitMolByPDBChainId",
    "SplitMolByPDBResidues",
    "StereoBondThresholds",
    "StereoDescriptor",
    "StereoGroup",
    "StereoGroupType",
    "StereoGroup_vect",
    "StereoInfo",
    "StereoSpecified",
    "StereoType",
    "SubstanceGroup",
    "SubstanceGroupAttach",
    "SubstanceGroupCState",
    "SubstanceGroup_VECT",
    "SubstructMatchParameters",
    "SupplierFromFilename",
    "TDTMolSupplier",
    "TDTWriter",
    "TranslateChiralFlagToStereoGroups",
    "UNCONSTRAINED_ANIONS",
    "UNCONSTRAINED_CATIONS",
    "UnfoldedRDKFingerprintCountBased",
    "WedgeBond",
    "WedgeMolBonds",
    "inchi",
    "molzip",
    "molzipFragments",
    "rdBase",
    "rdCIPLabeler",
    "rdCoordGen",
    "rdGeometry",
    "rdMolInterchange",
    "rdchem",
    "rdinchi",
    "rdmolfiles",
    "rdmolops",
    "templDir",
    "tossit"
]


class _GetRDKitObjIterator():
    pass
class _GetBondsIterator(_GetRDKitObjIterator):
    pass
class _GetAtomsIterator(_GetRDKitObjIterator):
    pass
def AddHs( mol: Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: object = None, addResidueInfo: bool = False) -> Mol:
    """
    AddHs( mol: Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: object = None, addResidueInfo: bool = False) -> Mol
        Adds hydrogens to the graph of a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - explicitOnly: (optional) if this toggle is set, only explicit Hs will
              be added to the molecule.  Default value is 0 (add implicit and explicit Hs).
        
            - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates
              set.  Default value is 0 (no 3D coords).
        
            - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
              considered to have Hs added to them
        
            - addResidueInfo: (optional) if this is true, add residue info to
              hydrogen atoms (useful for PDB files).
        
          RETURNS: a new molecule with added Hs
        
          NOTES:
        
            - The original molecule is *not* modified.
        
            - Much of the code assumes that Hs are not included in the molecular
              topology, so be *very* careful with the molecule that comes back from
              this function.
        
        

        C++ signature :
            RDKit::ROMol* AddHs(RDKit::ROMol [,bool=False [,bool=False [,boost::python::api::object=None [,bool=False]]]])
    """
def AddMetadataToPNGFile( metadata: dict, filename: object) -> object:
    """
    AddMetadataToPNGFile( metadata: dict, filename: object) -> object
        Adds metadata to PNG data read from a file.
        
             ARGUMENTS:
        
               - metadata: dict with the metadata to be written
                           (keys and values should be strings)
        
               - filename: the PNG filename
        
             RETURNS:
               the updated PNG data

        C++ signature :
            boost::python::api::object AddMetadataToPNGFile(boost::python::dict,boost::python::api::object)
    """
def AddMetadataToPNGString( metadata: dict, png: object) -> object:
    """
    AddMetadataToPNGString( metadata: dict, png: object) -> object
        Adds metadata to a PNG string.
        
             ARGUMENTS:
        
               - metadata: dict with the metadata to be written
                           (keys and values should be strings)
        
               - png: the PNG string
        
             RETURNS:
               the updated PNG data

        C++ signature :
            boost::python::api::object AddMetadataToPNGString(boost::python::dict,boost::python::api::object)
    """
def AddMolSubstanceGroup( mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup:
    """
    AddMolSubstanceGroup( mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup
        adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup

        C++ signature :
            RDKit::SubstanceGroup* AddMolSubstanceGroup(RDKit::ROMol {lvalue},RDKit::SubstanceGroup)
    """
def AddRecursiveQuery( mol: Mol, query: Mol, atomIdx: int, preserveExistingQuery: bool = True) -> None:
    """
    AddRecursiveQuery( mol: Mol, query: Mol, atomIdx: int, preserveExistingQuery: bool = True) -> None
        Adds a recursive query to an atom
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as the recursive query (this will be copied)
        
            - atomIdx: the atom to modify
        
            - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved
        
          RETURNS: None
        
        

        C++ signature :
            void AddRecursiveQuery(RDKit::ROMol {lvalue},RDKit::ROMol,unsigned int [,bool=True])
    """
def AddWavyBondsForStereoAny( mol: Mol, clearDoubleBondFlags: bool = True, addWhenImpossible: int = 1000) -> None:
    """
    AddWavyBondsForStereoAny( mol: Mol, clearDoubleBondFlags: bool = True, addWhenImpossible: int = 1000) -> None
        set wavy bonds around double bonds with STEREOANY stereo
          ARGUMENTS :
            - molecule : the molecule to update\n -
            - conformer : the conformer to use to determine wedge direction
        

        C++ signature :
            void AddWavyBondsForStereoAny(RDKit::ROMol {lvalue} [,bool=True [,unsigned int=1000]])
    """
def AdjustQueryProperties( mol: Mol, params: object = None) -> Mol:
    """
    AdjustQueryProperties( mol: Mol, params: object = None) -> Mol
        Returns a new molecule where the query properties of atoms have been modified.

        C++ signature :
            RDKit::ROMol* AdjustQueryProperties(RDKit::ROMol [,boost::python::api::object=None])
    """
def AdjustQueryPropertiesWithGenericGroups( mol: Mol, params: object = None) -> Mol:
    """
    AdjustQueryPropertiesWithGenericGroups( mol: Mol, params: object = None) -> Mol
        Returns a new molecule where the query properties of atoms have been modified and generic group queries have been prepared.

        C++ signature :
            RDKit::ROMol* AdjustQueryPropertiesWithGenericGroups(RDKit::ROMol [,boost::python::api::object=None])
    """
def AssignAtomChiralTagsFromMolParity( mol: Mol, replaceExistingTags: bool = True) -> None:
    """
    AssignAtomChiralTagsFromMolParity( mol: Mol, replaceExistingTags: bool = True) -> None
        Sets the chiral tags on a molecule's atoms based on
          the molParity atom property.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - replaceExistingTags: if True, existing stereochemistry information will be cleared
            before running the calculation.
        
        

        C++ signature :
            void AssignAtomChiralTagsFromMolParity(RDKit::ROMol {lvalue} [,bool=True])
    """
def AssignAtomChiralTagsFromStructure( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    AssignAtomChiralTagsFromStructure( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None
        Sets the chiral tags on a molecule's atoms based on
          a 3D conformation.
          NOTE that this does not check to see if atoms are chiral centers (i.e. all
          substituents are different), it merely sets the chiral type flags based on the
          coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you
          want chiral flags only on actual stereocenters.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: the conformer id to use, -1 for the default 
            - replaceExistingTags: if True, existing stereochemistry information will be cleared
            before running the calculation.
        
        

        C++ signature :
            void AssignAtomChiralTagsFromStructure(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
def AssignCIPLabels( mol: Mol, atomsToLabel: object = None, bondsToLabel: object = None, maxRecursiveIterations: int = 0) -> None:
    """
    AssignCIPLabels( mol: Mol, atomsToLabel: object = None, bondsToLabel: object = None, maxRecursiveIterations: int = 0) -> None
        New implementation of Stereo assignment using a true CIP ranking.
        On return:  The molecule to contains CIP flags
        Errors:  when maxRecursiveIterations is exceeded, throws a MaxIterationsExceeded error
        ARGUMENTS:
        
         - mol: the molecule
         - atomsToLabel: (optional) list of atoms to label
         - bondsToLabel: (optional) list of bonds to label
         - maxRecursiveIterations: (optional) protects against pseudo-infinite
        recursion for highly symmetrical structures.
         A value of 1,250,000 take about 1 second.  Most structures requires less than 10,000iterations.
         A peptide with MW~3000 took about 100 iterations, and a 20,000 mw protein took about 600 iterations
        (0 = default - no limit)
        

        C++ signature :
            void AssignCIPLabels(RDKit::ROMol {lvalue} [,boost::python::api::object=None [,boost::python::api::object=None [,unsigned int=0]]])
    """
def AssignChiralTypesFromBondDirs( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    AssignChiralTypesFromBondDirs( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None
        Uses bond directions to assign ChiralTypes to a molecule's atoms.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: (optional) the conformation to use 
            - replaceExistingTags: (optional) replace any existing information about stereochemistry
        
        

        C++ signature :
            void AssignChiralTypesFromBondDirs(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
def AssignRadicals( mol: Mol) -> None:
    """
    AssignRadicals( mol: Mol) -> None
        Assigns radical counts to atoms
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void AssignRadicals(RDKit::ROMol {lvalue})
    """
def AssignStereochemistry( mol: Mol, cleanIt: bool = False, force: bool = False, flagPossibleStereoCenters: bool = False) -> None:
    """
    AssignStereochemistry( mol: Mol, cleanIt: bool = False, force: bool = False, flagPossibleStereoCenters: bool = False) -> None
        Does the CIP stereochemistry assignment 
          for the molecule's atoms (R/S) and double bond (Z/E).
          Chiral atoms will have a property '_CIPCode' indicating
          their chiral code.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - cleanIt: (optional) if provided, any existing values of the property `_CIPCode`
                will be cleared, atoms with a chiral specifier that aren't
              actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
              etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with 
              STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP 
              atom ranks will be marked STEREONONE. 
            - force: (optional) causes the calculation to be repeated, even if it has already
              been done
            - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
              atoms that are possible stereocenters
        

        C++ signature :
            void AssignStereochemistry(RDKit::ROMol {lvalue} [,bool=False [,bool=False [,bool=False]]])
    """
def AssignStereochemistryFrom3D( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
    AssignStereochemistryFrom3D( mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None
        Uses a conformer (should be 3D) to assign ChiralTypes to a molecule's atoms
                and stereo flags to its bonds
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: (optional) the conformation to use 
            - replaceExistingTags: (optional) replace any existing information about stereochemistry
        
        

        C++ signature :
            void AssignStereochemistryFrom3D(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
def AtomFromSmarts( SMARTS: str) -> Atom:
    """
    AtomFromSmarts( SMARTS: str) -> Atom
        Construct an atom from a SMARTS string

        C++ signature :
            RDKit::Atom* AtomFromSmarts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def AtomFromSmiles( SMILES: str) -> Atom:
    """
    AtomFromSmiles( SMILES: str) -> Atom
        Construct an atom from a SMILES string

        C++ signature :
            RDKit::Atom* AtomFromSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def BondFromSmarts( SMILES: str) -> Bond:
    """
    BondFromSmarts( SMILES: str) -> Bond
        Construct a bond from a SMARTS string

        C++ signature :
            RDKit::Bond* BondFromSmarts(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def BondFromSmiles( SMILES: str) -> Bond:
    """
    BondFromSmiles( SMILES: str) -> Bond
        Construct a bond from a SMILES string

        C++ signature :
            RDKit::Bond* BondFromSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def CanonicalRankAtoms( mol: Mol, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vectj:
    """
    CanonicalRankAtoms( mol: Mol, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vectj
        Returns the canonical atom ranking for each atom of a molecule fragment.
          If breakTies is False, this returns the symmetry class for each atom.  The symmetry
          class is used by the canonicalization routines to type each atom based on the whole
          chemistry of the molecular graph.  Any atom with the same rank (symmetry class) is
          indistinguishable.  For example:
        
            >>> mol = MolFromSmiles('C1NCN1')
            >>> list(CanonicalRankAtoms(mol, breakTies=False))
            [0,1,0,1]
        
          In this case the carbons have the same symmetry class and the nitrogens have the same
          symmetry class.  From the perspective of the Molecular Graph, they are identical.
        
          ARGUMENTS:
        
            - mol: the molecule
            - breakTies: (optional) force breaking of ranked ties [default=True]
            - includeChirality: (optional) use chiral information when computing rank [default=True]
            - includeIsotopes: (optional) use isotope information when computing rank [default=True]
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::vector<unsigned int, std::allocator<unsigned int> > CanonicalRankAtoms(RDKit::ROMol [,bool=True [,bool=True [,bool=True]]])
    """
def CanonicalRankAtomsInFragment( mol: Mol, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vecti:
    """
    CanonicalRankAtomsInFragment( mol: Mol, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vecti
        Returns the canonical atom ranking for each atom of a molecule fragment
          See help(CanonicalRankAtoms) for more information.
        
           >>> mol = MolFromSmiles('C1NCN1.C1NCN1')
           >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(0,4), breakTies=False))
           [4,6,4,6,-1,-1,-1,-1]
           >>> list(CanonicalRankAtomsInFragment(mol, atomsToUse=range(4,8), breakTies=False))
           [-1,-1,-1,-1,4,6,4,6]
        
          ARGUMENTS:
        
            - mol: the molecule
            - atomsToUse : a list of atoms to include in the fragment
            - bondsToUse : (optional) a list of bonds to include in the fragment
              if not provided, all bonds between the atoms provided
              will be included.
            - atomSymbols : (optional) a list with the symbols to use for the atoms
              in the SMILES. This should have be mol.GetNumAtoms() long.
            - breakTies: (optional) force breaking of ranked ties
            - includeChirality: (optional) use chiral information when computing rank [default=True]
            - includeIsotopes: (optional) use isotope information when computing rank [default=True]
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::vector<int, std::allocator<int> > CanonicalRankAtomsInFragment(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=True [,bool=True]]]]])
    """
def CanonicalizeEnhancedStereo( mol: Mol) -> None:
    """
    CanonicalizeEnhancedStereo( mol: Mol) -> None

        C++ signature :
            void CanonicalizeEnhancedStereo(RDKit::ROMol {lvalue})
    """
def Cleanup( mol: Mol) -> None:
    """
    Cleanup( mol: Mol) -> None
        cleans up certain common bad functionalities in the molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void Cleanup(RDKit::ROMol {lvalue})
    """
def CleanupOrganometallics( mol: Mol) -> None:
    """
    CleanupOrganometallics( mol: Mol) -> None
        cleans up certain common bad functionalities in the organometallic molecule
        
          Note that this function is experimental and may either change in behavior
          or be replaced with something else in future releases.
        
                ARGUMENTS :
        
         - mol : the molecule to use
        
         NOTES :
        
         - The molecule is modified in place.
        
         

        C++ signature :
            void CleanupOrganometallics(RDKit::ROMol {lvalue})
    """
def ClearMolSubstanceGroups( arg1: Mol) -> None:
    """
    ClearMolSubstanceGroups( arg1: Mol) -> None
        removes all SubstanceGroups from a molecule (if any)

        C++ signature :
            void ClearMolSubstanceGroups(RDKit::ROMol {lvalue})
    """
def CombineMols( mol1: Mol, mol2: Mol, offset: Point3D = Point3D()) -> Mol:
    """
    CombineMols( mol1: Mol, mol2: Mol, offset: Point3D = Point3D()) -> Mol
        Combine the atoms from two molecules to produce a third

        C++ signature :
            RDKit::ROMol* CombineMols(RDKit::ROMol,RDKit::ROMol [,RDGeom::Point3D=<rdkit.Geometry.rdGeometry.Point3D object at 0x7f2c653ba1c0>])
    """
def ConvertGenericQueriesToSubstanceGroups( mol: Mol) -> None:
    """
    ConvertGenericQueriesToSubstanceGroups( mol: Mol) -> None
        documentation

        C++ signature :
            void ConvertGenericQueriesToSubstanceGroups(RDKit::ROMol {lvalue})
    """
def CreateAtomBoolPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    CreateAtomBoolPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None
        creates a list property on the molecule from individual atom property values

        C++ signature :
            void CreateAtomBoolPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,unsigned int=190]])
    """
def CreateAtomDoublePropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    CreateAtomDoublePropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None
        creates a list property on the molecule from individual atom property values

        C++ signature :
            void CreateAtomDoublePropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,unsigned int=190]])
    """
def CreateAtomIntPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    CreateAtomIntPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None
        creates a list property on the molecule from individual atom property values

        C++ signature :
            void CreateAtomIntPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,unsigned int=190]])
    """
def CreateAtomStringPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None:
    """
    CreateAtomStringPropertyList( mol: Mol, propName: str, missingValueMarker: str = '', lineSize: int = 190) -> None
        creates a list property on the molecule from individual atom property values

        C++ signature :
            void CreateAtomStringPropertyList(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,unsigned int=190]])
    """
def CreateMolDataSubstanceGroup( mol: Mol, fieldName: str, value: str) -> SubstanceGroup:
    """
    CreateMolDataSubstanceGroup( mol: Mol, fieldName: str, value: str) -> SubstanceGroup
        creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup

        C++ signature :
            RDKit::SubstanceGroup* CreateMolDataSubstanceGroup(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def CreateMolSubstanceGroup( mol: Mol, type: str) -> SubstanceGroup:
    """
    CreateMolSubstanceGroup( mol: Mol, type: str) -> SubstanceGroup
        creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup

        C++ signature :
            RDKit::SubstanceGroup* CreateMolSubstanceGroup(RDKit::ROMol {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def CreateStereoGroup( stereoGroupType: StereoGroupType, mol: Mol, atomIds: object, readId: int = 0) -> StereoGroup:
    """
    CreateStereoGroup( stereoGroupType: StereoGroupType, mol: Mol, atomIds: object, readId: int = 0) -> StereoGroup
        creates a StereoGroup associated with a molecule from a list of atom Ids

        C++ signature :
            RDKit::StereoGroup* CreateStereoGroup(RDKit::StereoGroupType,RDKit::ROMol {lvalue},boost::python::api::object [,unsigned int=0])
    """
def DativeBondsToHaptic( mol: Mol) -> Mol:
    """
    DativeBondsToHaptic( mol: Mol) -> Mol
        Does the reverse of hapticBondsToDative.  If there are multiple
        contiguous atoms attached by dative bonds to an atom (probably a metal
        atom), the dative bonds will be replaced by a dummy atom in their
        centre attached to the (metal) atom by a dative bond, which is
        labelled with ENDPTS of the atoms that had the original dative bonds.
        
        ARGUMENTS:
        
          - mol: the molecule to use
        
        RETURNS:
          a modified copy of the molecule

        C++ signature :
            RDKit::ROMol* DativeBondsToHaptic(RDKit::ROMol)
    """
def DeleteSubstructs( mol: Mol, query: Mol, onlyFrags: bool = False, useChirality: bool = False) -> Mol:
    """
    DeleteSubstructs( mol: Mol, query: Mol, onlyFrags: bool = False, useChirality: bool = False) -> Mol
        Removes atoms matching a substructure query from a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as a substructure query
        
            - onlyFrags: (optional) if this toggle is set, atoms will only be removed if
              the entire fragment in which they are found is matched by the query.
              See below for examples.
              Default value is 0 (remove the atoms whether or not the entire fragment matches)
        
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a new molecule with the substructure removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - DeleteSubstructs('CCOC','OC') -> 'CC'
        
            - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'
        
            - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'
        
            - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'
        
        

        C++ signature :
            RDKit::ROMol* DeleteSubstructs(RDKit::ROMol,RDKit::ROMol [,bool=False [,bool=False]])
    """
def DetectBondStereoChemistry( mol: Mol, conformer: Conformer) -> None:
    """
    DetectBondStereoChemistry( mol: Mol, conformer: Conformer) -> None
        Assign stereochemistry to bonds based on coordinates and a conformer.
                DEPRECATED
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
            - conformer: Conformer providing the coordinates
        
        

        C++ signature :
            void DetectBondStereoChemistry(RDKit::ROMol {lvalue},RDKit::Conformer const*)
    """
def DetectBondStereochemistry( mol: Mol, confId: int = -1) -> None:
    """
    DetectBondStereochemistry( mol: Mol, confId: int = -1) -> None
        DEPRECATED, use SetDoubleBondNeighborDirections() instead
          ARGUMENTS:
          
            - mol: the molecule to be modified
            - confId: Conformer to use for the coordinates
        
        

        C++ signature :
            void DetectBondStereochemistry(RDKit::ROMol {lvalue} [,int=-1])
    """
def DetectChemistryProblems( mol: Mol, sanitizeOps: int = rdmolops.SanitizeFlags.SANITIZE_ALL) -> tuple:
    """
    DetectChemistryProblems( mol: Mol, sanitizeOps: int = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL) -> tuple
        checks for chemistry problems

        C++ signature :
            boost::python::tuple DetectChemistryProblems(RDKit::ROMol [,unsigned int=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL])
    """
def FastFindRings( arg1: Mol) -> None:
    """
    FastFindRings( arg1: Mol) -> None
        Does a non-SSSR ring finding for a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use.
        
          RETURNS: Nothing
        
        

        C++ signature :
            void FastFindRings(RDKit::ROMol)
    """
def FindAllPathsOfLengthN( mol: Mol, length: int, useBonds: bool = True, useHs: bool = False, rootedAtAtom: int = -1, onlyShortestPaths: bool = False) -> _listSt6vectorIiSaIiEE:
    """
    FindAllPathsOfLengthN( mol: Mol, length: int, useBonds: bool = True, useHs: bool = False, rootedAtAtom: int = -1, onlyShortestPaths: bool = False) -> _listSt6vectorIiSaIiEE
        Finds all paths of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target length for the paths.
        
            - useBonds: (optional) toggles the use of bond indices in the paths.
              Otherwise atom indices are used.  *Note* this behavior is different
              from that for subgraphs.
              Defaults to 1.
        
            - rootedAtAtom: (optional) if nonzero, only paths from the specified
              atom will be returned.
        
            - onlyShortestPaths: (optional) if set then only paths which are <= the shortest
              path between the begin and end atoms will be included in the results
        
          RETURNS: a tuple of tuples with IDs for the bonds.
        
          NOTES: 
        
           - Difference between _subgraphs_ and _paths_ :: 
        
               Subgraphs are potentially branched, whereas paths (in our 
               terminology at least) cannot be.  So, the following graph: 
        
                    C--0--C--1--C--3--C
                          |
                          2
                          |
                          C
        
               has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
               but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
        
        

        C++ signature :
            std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindAllPathsOfLengthN(RDKit::ROMol,unsigned int [,bool=True [,bool=False [,int=-1 [,bool=False]]]])
    """
def FindAllSubgraphsOfLengthMToN( mol: Mol, min: int, max: int, useHs: bool = False, rootedAtAtom: int = -1) -> object:
    """
    FindAllSubgraphsOfLengthMToN( mol: Mol, min: int, max: int, useHs: bool = False, rootedAtAtom: int = -1) -> object
        Finds all subgraphs of a particular length in a molecule
          See documentation for FindAllSubgraphsOfLengthN for definitions
        
        

        C++ signature :
            boost::python::api::object FindAllSubgraphsOfLengthMToN(RDKit::ROMol,unsigned int,unsigned int [,bool=False [,int=-1]])
    """
def FindAllSubgraphsOfLengthN( mol: Mol, length: int, useHs: bool = False, rootedAtAtom: int = -1) -> _listSt6vectorIiSaIiEE:
    """
    FindAllSubgraphsOfLengthN( mol: Mol, length: int, useHs: bool = False, rootedAtAtom: int = -1) -> _listSt6vectorIiSaIiEE
        Finds all subgraphs of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target number of bonds for the subgraphs.
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
              atom will be returned.
        
          RETURNS: a tuple of 2-tuples with bond IDs
        
          NOTES: 
        
           - Difference between _subgraphs_ and _paths_ :: 
        
               Subgraphs are potentially branched, whereas paths (in our 
               terminology at least) cannot be.  So, the following graph: 
        
                    C--0--C--1--C--3--C
                          |
                          2
                          |
                          C
          has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
          but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
        
        

        C++ signature :
            std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindAllSubgraphsOfLengthN(RDKit::ROMol,unsigned int [,bool=False [,int=-1]])
    """
def FindAtomEnvironmentOfRadiusN( mol: Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: object = None) -> _vecti:
    """
    FindAtomEnvironmentOfRadiusN( mol: Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: object = None) -> _vecti
        Find bonds of a particular radius around an atom. 
                 Return empty result if there is no bond at the requested radius.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - radius: an integer with the target radius for the environment.
        
            - rootedAtAtom: the atom to consider
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - enforceSize (optional) If set to False, all bonds within the requested radius is 
              collected. Defaults to 1. 
        
            - atomMap: (optional) If provided, it will measure the minimum distance of the atom 
              from the rooted atom (start with 0 from the rooted atom). The result is a pair of 
              the atom ID and the distance. 
        
          RETURNS: a vector of bond IDs
        
        

        C++ signature :
            std::vector<int, std::allocator<int> > FindAtomEnvironmentOfRadiusN(RDKit::ROMol,unsigned int,unsigned int [,bool=False [,bool=True [,boost::python::api::object=None]]])
    """
def FindPotentialStereo( mol: Mol, cleanIt: bool = False, flagPossible: bool = True) -> _vectN5RDKit9Chirality10StereoInfoE:
    """
    FindPotentialStereo( mol: Mol, cleanIt: bool = False, flagPossible: bool = True) -> _vectN5RDKit9Chirality10StereoInfoE
        find potential stereo elements in a molecule and returns them as StereoInfo objects
        Note that this function is still somewhat experimental and the API
        and results may change in a future release.

        C++ signature :
            std::vector<RDKit::Chirality::StereoInfo, std::allocator<RDKit::Chirality::StereoInfo> > FindPotentialStereo(RDKit::ROMol {lvalue} [,bool=False [,bool=True]])
    """
def FindPotentialStereoBonds( mol: Mol, cleanIt: bool = False) -> None:
    """
    FindPotentialStereoBonds( mol: Mol, cleanIt: bool = False) -> None
        Find bonds than can be cis/trans in a molecule and mark them as 'any'.
                 This function finds any double bonds that can potentially be part
                 of a cis/trans system. No attempt is made here to mark them cis or trans
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode
                       on the bond is cleared - otherwise it is left untouched
        
        

        C++ signature :
            void FindPotentialStereoBonds(RDKit::ROMol {lvalue} [,bool=False])
    """
def FindRingFamilies( arg1: Mol) -> None:
    """
    FindRingFamilies( arg1: Mol) -> None
        generate Unique Ring Families

        C++ signature :
            void FindRingFamilies(RDKit::ROMol)
    """
def FindUniqueSubgraphsOfLengthN( mol: Mol, length: int, useHs: bool = False, useBO: bool = True, rootedAtAtom: int = -1) -> _listSt6vectorIiSaIiEE:
    """
    FindUniqueSubgraphsOfLengthN( mol: Mol, length: int, useHs: bool = False, useBO: bool = True, rootedAtAtom: int = -1) -> _listSt6vectorIiSaIiEE
        Finds unique subgraphs of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target number of bonds for the subgraphs.
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from
              another.
              Defaults to 1.
        
            - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
              atom will be returned.
        
          RETURNS: a tuple of tuples with bond IDs
        
        
        

        C++ signature :
            std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindUniqueSubgraphsOfLengthN(RDKit::ROMol,unsigned int [,bool=False [,bool=True [,int=-1]]])
    """
def ForwardStereoGroupIds( arg1: Mol) -> None:
    """
    ForwardStereoGroupIds( arg1: Mol) -> None
        Forward the original Stereo Group IDs when exporting the Mol.

        C++ signature :
            void ForwardStereoGroupIds(RDKit::ROMol {lvalue})
    """
def FragmentOnBRICSBonds( mol: Mol) -> Mol:
    """
    FragmentOnBRICSBonds( mol: Mol) -> Mol
        Return a new molecule with all BRICS bonds broken

        C++ signature :
            RDKit::ROMol* FragmentOnBRICSBonds(RDKit::ROMol)
    """
def FragmentOnBonds( mol: Mol, bondIndices: object, addDummies: bool = True, dummyLabels: object = None, bondTypes: object = None, cutsPerAtom: list = []) -> Mol:
    """
    FragmentOnBonds( mol: Mol, bondIndices: object, addDummies: bool = True, dummyLabels: object = None, bondTypes: object = None, cutsPerAtom: list = []) -> Mol
        Return a new molecule with all specified bonds broken
        
          ARGUMENTS:
        
              - mol            - the molecule to be modified
              - bondIndices    - indices of the bonds to be broken
              - addDummies  - toggles addition of dummy atoms to indicate where 
                bonds were broken
              - dummyLabels - used to provide the labels to be used for the dummies.
                the first element in each pair is the label for the dummy
                that replaces the bond's beginAtom, the second is for the 
                dummy that replaces the bond's endAtom. If not provided, the
                dummies are labeled with atom indices.
              - bondTypes - used to provide the bond type to use between the
                fragments and the dummy atoms. If not provided, defaults to single. 
              - cutsPerAtom - used to return the number of cuts made at each atom. 
        
          RETURNS:
              a new Mol with the modifications
        

        C++ signature :
            RDKit::ROMol* FragmentOnBonds(RDKit::ROMol,boost::python::api::object [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::list=[]]]]])
    """
def FragmentOnSomeBonds( mol: Mol, bondIndices: object, numToBreak: int = 1, addDummies: bool = True, dummyLabels: object = None, bondTypes: object = None, returnCutsPerAtom: bool = False) -> tuple:
    """
    FragmentOnSomeBonds( mol: Mol, bondIndices: object, numToBreak: int = 1, addDummies: bool = True, dummyLabels: object = None, bondTypes: object = None, returnCutsPerAtom: bool = False) -> tuple
        fragment on some bonds

        C++ signature :
            boost::python::tuple FragmentOnSomeBonds(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None [,bool=False]]]]])
    """
def Get3DDistanceMatrix( mol: Mol, confId: int = -1, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> object:
    """
    Get3DDistanceMatrix( mol: Mol, confId: int = -1, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> object
        Returns the molecule's 3D distance matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - confId: (optional) chooses the conformer Id to use
              Default value is -1.
        
            - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
              matrix (to return a "Balaban" distance matrix).
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats with the distance matrix
        
        

        C++ signature :
            _object* Get3DDistanceMatrix(RDKit::ROMol {lvalue} [,int=-1 [,bool=False [,bool=False [,char const*='']]]])
    """
def GetAdjacencyMatrix( mol: Mol, useBO: bool = False, emptyVal: int = 0, force: bool = False, prefix: str = '') -> object:
    """
    GetAdjacencyMatrix( mol: Mol, useBO: bool = False, emptyVal: int = 0, force: bool = False, prefix: str = '') -> object
        Returns the molecule's adjacency matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - useBO: (optional) toggles use of bond orders in calculating the matrix.
              Default value is 0.
        
            - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats containing the adjacency matrix
        
        

        C++ signature :
            _object* GetAdjacencyMatrix(RDKit::ROMol {lvalue} [,bool=False [,int=0 [,bool=False [,char const*='']]]])
    """
def GetAllowNontetrahedralChirality() -> bool:
    """
    GetAllowNontetrahedralChirality() -> bool
        returns whether or not recognition of non-tetrahedral chirality from 3D structures is enabled

        C++ signature :
            bool GetAllowNontetrahedralChirality()
    """
def GetAtomAlias( atom: Atom) -> str:
    """
    GetAtomAlias( atom: Atom) -> str
        Returns the atom's MDL alias text

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAtomAlias(RDKit::Atom const*)
    """
def GetAtomRLabel( atom: Atom) -> int:
    """
    GetAtomRLabel( atom: Atom) -> int
        Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)

        C++ signature :
            int GetAtomRLabel(RDKit::Atom const*)
    """
def GetAtomValue( atom: Atom) -> str:
    """
    GetAtomValue( atom: Atom) -> str
        Returns the atom's MDL alias text

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAtomValue(RDKit::Atom const*)
    """
def GetDefaultPickleProperties() -> int:
    """
    GetDefaultPickleProperties() -> int
        Get the current global mol pickler options.

        C++ signature :
            unsigned int GetDefaultPickleProperties()
    """
def GetDistanceMatrix( mol: Mol, useBO: bool = False, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> object:
    """
    GetDistanceMatrix( mol: Mol, useBO: bool = False, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> object
        Returns the molecule's topological distance matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - useBO: (optional) toggles use of bond orders in calculating the distance matrix.
              Default value is 0.
        
            - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
              matrix (to return a "Balaban" distance matrix).
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats with the distance matrix
        
        

        C++ signature :
            _object* GetDistanceMatrix(RDKit::ROMol {lvalue} [,bool=False [,bool=False [,bool=False [,char const*='']]]])
    """
def GetFormalCharge( arg1: Mol) -> int:
    """
    GetFormalCharge( arg1: Mol) -> int
        Returns the formal charge for the molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
        

        C++ signature :
            int GetFormalCharge(RDKit::ROMol)
    """
def GetMolFrags( mol: Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: object = None, fragsMolAtomMapping: object = None) -> tuple:
    """
    GetMolFrags( mol: Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: object = None, fragsMolAtomMapping: object = None) -> tuple
        Finds the disconnected fragments from a molecule.
        
          For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns
          ((0, 1, 2, 3), (4, 5))
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - asMols: (optional) if this is provided and true, the fragments
              will be returned as molecules instead of atom ids.
            - sanitizeFrags: (optional) if this is provided and true, the fragments
              molecules will be sanitized before returning them.
            - frags: (optional, defaults to None) if asMols is true and this is provided
               as an empty list, the result will be mol.GetNumAtoms() long on return and
               will contain the fragment assignment for each Atom
            - fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this
              is provided as an empty list, the result will be numFrags long on 
              return, and each entry will contain the indices of the Atoms in that fragment:
              [(0, 1, 2, 3), (4, 5)]
        
          RETURNS: a tuple of tuples with IDs for the atoms in each fragment
                   or a tuple of molecules.
        
        

        C++ signature :
            boost::python::tuple GetMolFrags(RDKit::ROMol [,bool=False [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None]]]])
    """
def GetMolSubstanceGroupWithIdx( arg1: Mol, arg2: int) -> SubstanceGroup:
    """
    GetMolSubstanceGroupWithIdx( arg1: Mol, arg2: int) -> SubstanceGroup
        returns a particular SubstanceGroup from the molecule

        C++ signature :
            RDKit::SubstanceGroup* GetMolSubstanceGroupWithIdx(RDKit::ROMol {lvalue},unsigned int)
    """
def GetMolSubstanceGroups( arg1: Mol) -> SubstanceGroup_VECT:
    """
    GetMolSubstanceGroups( arg1: Mol) -> SubstanceGroup_VECT
        returns a copy of the molecule's SubstanceGroups (if any)

        C++ signature :
            std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > GetMolSubstanceGroups(RDKit::ROMol {lvalue})
    """
def GetMostSubstitutedCoreMatch( mol: Mol, core: Mol, matches: object) -> object:
    """
    GetMostSubstitutedCoreMatch( mol: Mol, core: Mol, matches: object) -> object
        Postprocesses the results of a mol.GetSubstructMatches(core) call 
        where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
        It returns the match with the largest number of non-hydrogen matches to 
        the terminal dummy atoms.
        
          ARGUMENTS:
        
            - mol: the molecule GetSubstructMatches was run on
        
            - core: the molecule used as a substructure query
        
            - matches: the result returned by GetSubstructMatches
        
          RETURNS: the tuple where terminal dummy atoms in the core match the largest 
                   number of non-hydrogen atoms in mol
        

        C++ signature :
            _object* GetMostSubstitutedCoreMatch(RDKit::ROMol,RDKit::ROMol,boost::python::api::object)
    """
def GetPeriodicTable() -> PeriodicTable:
    """
    GetPeriodicTable() -> PeriodicTable
        Returns the application's PeriodicTable instance.
        
        

        C++ signature :
            RDKit::PeriodicTable* GetPeriodicTable()
    """
def GetSSSR( mol: Mol, includeDativeBonds: bool = False) -> _vectSt6vectorIiSaIiEE:
    """
    GetSSSR( mol: Mol, includeDativeBonds: bool = False) -> _vectSt6vectorIiSaIiEE
        Get the smallest set of simple rings for a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use.
            - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
        
          RETURNS: a sequence of sequences containing the rings found as atom ids
                 The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.
        
        

        C++ signature :
            std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > GetSSSR(RDKit::ROMol {lvalue} [,bool=False])
    """
def GetShortestPath( arg1: Mol, arg2: int, arg3: int) -> tuple:
    """
    GetShortestPath( arg1: Mol, arg2: int, arg3: int) -> tuple
        Find the shortest path between two atoms using the Bellman-Ford algorithm.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - idx1: index of the first atom
            - idx2: index of the second atom
        
        

        C++ signature :
            boost::python::tuple GetShortestPath(RDKit::ROMol,int,int)
    """
def GetSupplementalSmilesLabel( atom: Atom) -> str:
    """
    GetSupplementalSmilesLabel( atom: Atom) -> str
        Gets the supplemental smiles label on an atom, returns an empty string if not present.

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSupplementalSmilesLabel(RDKit::Atom const*)
    """
def GetSymmSSSR( mol: Mol, includeDativeBonds: bool = False) -> _vectSt6vectorIiSaIiEE:
    """
    GetSymmSSSR( mol: Mol, includeDativeBonds: bool = False) -> _vectSt6vectorIiSaIiEE
        Get a symmetrized SSSR for a molecule.
        
          The symmetrized SSSR is at least as large as the SSSR for a molecule.
          In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be
          a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).
        
          ARGUMENTS:
        
            - mol: the molecule to use.
            - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
        
          RETURNS: a sequence of sequences containing the rings found as atom ids
        
        

        C++ signature :
            std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > GetSymmSSSR(RDKit::ROMol {lvalue} [,bool=False])
    """
def GetUseLegacyStereoPerception() -> bool:
    """
    GetUseLegacyStereoPerception() -> bool
        returns whether or not the legacy stereo perception code is being used

        C++ signature :
            bool GetUseLegacyStereoPerception()
    """
def HapticBondsToDative( mol: Mol) -> Mol:
    """
    HapticBondsToDative( mol: Mol) -> Mol
        One way of showing haptic bonds (such as cyclopentadiene to
        iron in ferrocene) is to use a dummy atom with a dative bond to the
        iron atom with the bond labelled with the atoms involved in the
        organic end of the bond.  Another way is to have explicit dative
        bonds from the atoms of the haptic group to the metal atom.  This
        function converts the former representation to the latter.
        
        ARGUMENTS:
        
          - mol: the molecule to use
        
        RETURNS:
          a modified copy of the molecule

        C++ signature :
            RDKit::ROMol* HapticBondsToDative(RDKit::ROMol)
    """
def JSONToMols( jsonBlock: str, params: object = None) -> tuple:
    """
    JSONToMols( jsonBlock: str, params: object = None) -> tuple
        Convert JSON to a tuple of molecules
        
            ARGUMENTS:
              - jsonBlock: the molecule to work with
              - params: (optional) JSONParseParameters controlling the JSON parsing
            RETURNS:
              a tuple of Mols
        

        C++ signature :
            boost::python::tuple JSONToMols(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=None])
    """
def Kekulize( mol: Mol, clearAromaticFlags: bool = False) -> None:
    """
    Kekulize( mol: Mol, clearAromaticFlags: bool = False) -> None
        Kekulizes the molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
              molecule will be marked non-aromatic following the kekulization.
              Default value is False.
        
          NOTES:
        
            - The molecule is modified in place.
        
            - this does not modify query bonds which have bond type queries (like those
              which come from SMARTS) or rings containing them.
        
            - even if clearAromaticFlags is False the BondType for all modified
              aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
              Kekulization.
        
        

        C++ signature :
            void Kekulize(RDKit::ROMol {lvalue} [,bool=False])
    """
def KekulizeIfPossible( mol: Mol, clearAromaticFlags: bool = False) -> None:
    """
    KekulizeIfPossible( mol: Mol, clearAromaticFlags: bool = False) -> None
        Kekulizes the molecule if possible. Otherwise the molecule is not modified
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the 
              molecule will be marked non-aromatic if the kekulization succeds.
              Default value is False.
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void KekulizeIfPossible(RDKit::ROMol {lvalue} [,bool=False])
    """
def LayeredFingerprint( mol: Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, branchedPaths: bool = True, fromAtoms: object = 0) -> ExplicitBitVect:
    """
    LayeredFingerprint( mol: Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, branchedPaths: bool = True, fromAtoms: object = 0) -> ExplicitBitVect
        Returns a layered fingerprint for a molecule
        
          NOTE: This function is experimental. The API or results may change from
            release to release.
        
          Explanation of the algorithm below.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - layerFlags: (optional) which layers to include in the fingerprint
              See below for definitions. Defaults to all.
        
            - minPath: (optional) minimum number of bonds to include in the subgraphs
              Defaults to 1.
        
            - maxPath: (optional) maximum number of bonds to include in the subgraphs
              Defaults to 7.
        
            - fpSize: (optional) number of bits in the fingerprint
              Defaults to 2048.
        
            - atomCounts: (optional) 
              if provided, this should be a list at least as long as the number of atoms
              in the molecule. It will be used to provide the count of the number 
              of paths that set bits each atom is involved in.
              NOTE: the list is not zeroed out here.
        
            - setOnlyBits: (optional) 
              if provided, only bits that are set in this bit vector will be set
              in the result. This is essentially the same as doing:
                   res &= setOnlyBits
              but also has an impact on the atomCounts (if being used)
        
            - branchedPaths: (optional) if set both branched and unbranched paths will be
              used in the fingerprint.
              Defaults to True.
        
            - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
              starting from these atoms will be used.
              Defaults to empty.
        
          RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits
        
          Layer definitions:
             - 0x01: pure topology
             - 0x02: bond order
             - 0x04: atom types
             - 0x08: presence of rings
             - 0x10: ring sizes
             - 0x20: aromaticity
        
        
        

        C++ signature :
            ExplicitBitVect* LayeredFingerprint(RDKit::ROMol [,unsigned int=4294967295 [,unsigned int=1 [,unsigned int=7 [,unsigned int=2048 [,boost::python::list=[] [,ExplicitBitVect*=None [,bool=True [,boost::python::api::object=0]]]]]]]])
    """
def MergeQueryHs( mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> Mol:
    """
    MergeQueryHs( mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> Mol
        merges hydrogens into their neighboring atoms as queries

        C++ signature :
            RDKit::ROMol* MergeQueryHs(RDKit::ROMol [,bool=False [,bool=False]])
    """
def MetadataFromPNGFile( filename: object) -> dict:
    """
    MetadataFromPNGFile( filename: object) -> dict
        Returns a dict with all metadata from the PNG file. Keys are strings, values are bytes.

        C++ signature :
            boost::python::dict MetadataFromPNGFile(boost::python::api::object)
    """
def MetadataFromPNGString( png: object) -> dict:
    """
    MetadataFromPNGString( png: object) -> dict
        Returns a dict with all metadata from the PNG string. Keys are strings, values are bytes.

        C++ signature :
            boost::python::dict MetadataFromPNGString(boost::python::api::object)
    """
def MolAddRecursiveQueries( mol: Mol, queries: dict, propName: str) -> None:
    """
    MolAddRecursiveQueries( mol: Mol, queries: dict, propName: str) -> None
        Adds named recursive queries to atoms
        

        C++ signature :
            void MolAddRecursiveQueries(RDKit::ROMol {lvalue},boost::python::dict,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def MolBundleCanSerialize() -> bool:
    """
    MolBundleCanSerialize() -> bool
        Returns True if the MolBundle is serializable (requires boost serialization

        C++ signature :
            bool MolBundleCanSerialize()
    """
def MolFragmentToCXSmarts( mol: Mol, atomsToUse: object, bondsToUse: object = 0, isomericSmarts: bool = True) -> str:
    """
    MolFragmentToCXSmarts( mol: Mol, atomsToUse: object, bondsToUse: object = 0, isomericSmarts: bool = True) -> str
        Returns a SMARTS string for a fragment of a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - atomsToUse: indices of atoms to include in the SMARTS string
            - bondsToUse: indices of bonds to include in the SMARTS string (optional)
            - isomericSmarts: (optional) include information about stereochemistry in
              the SMARTS.  Defaults to true.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToCXSmarts(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,bool=True]])
    """
@typing.overload
def MolFragmentToCXSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0) -> str:
    """
    MolFragmentToCXSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0) -> str
        Returns the CXSMILES string for a fragment of a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - params: the SmilesWriteParams 
            - atomsToUse : a list of atoms to include in the fragment
            - bondsToUse : (optional) a list of bonds to include in the fragment
              if not provided, all bonds between the atoms provided
              will be included.
            - atomSymbols : (optional) a list with the symbols to use for the atoms
              in the SMILES. This should have be mol.GetNumAtoms() long.
            - bondSymbols : (optional) a list with the symbols to use for the bonds
              in the SMILES. This should have be mol.GetNumBonds() long.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToCXSmiles(RDKit::ROMol,RDKit::SmilesWriteParams,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0]]])

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToCXSmiles(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False]]]]]]]]])
    """
@typing.overload
def MolFragmentToCXSmiles( mol: Mol, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    pass
def MolFragmentToSmarts( mol: Mol, atomsToUse: object, bondsToUse: object = 0, isomericSmarts: bool = True) -> str:
    """
    MolFragmentToSmarts( mol: Mol, atomsToUse: object, bondsToUse: object = 0, isomericSmarts: bool = True) -> str
        Returns a SMARTS string for a fragment of a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - atomsToUse: indices of atoms to include in the SMARTS string
            - bondsToUse: indices of bonds to include in the SMARTS string (optional)
            - isomericSmarts: (optional) include information about stereochemistry in
              the SMARTS.  Defaults to true.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToSmarts(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,bool=True]])
    """
@typing.overload
def MolFragmentToSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0) -> str:
    """
    MolFragmentToSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0) -> str
        Returns the canonical SMILES string for a fragment of a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - params: the SmilesWriteParams 
            - atomsToUse : a list of atoms to include in the fragment
            - bondsToUse : (optional) a list of bonds to include in the fragment
              if not provided, all bonds between the atoms provided
              will be included.
            - atomSymbols : (optional) a list with the symbols to use for the atoms
              in the SMILES. This should have be mol.GetNumAtoms() long.
            - bondSymbols : (optional) a list with the symbols to use for the bonds
              in the SMILES. This should have be mol.GetNumBonds() long.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToSmiles(RDKit::ROMol,RDKit::SmilesWriteParams,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0]]])

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolFragmentToSmiles(RDKit::ROMol,boost::python::api::object [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False]]]]]]]]])
    """
@typing.overload
def MolFragmentToSmiles( mol: Mol, atomsToUse: object, bondsToUse: object = 0, atomSymbols: object = 0, bondSymbols: object = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    pass
def MolFromFASTA( text: object, sanitize: bool = True, flavor: int = 0) -> Mol:
    """
    MolFromFASTA( text: object, sanitize: bool = True, flavor: int = 0) -> Mol
        Construct a molecule from a FASTA string (currently only supports peptides).
        
          ARGUMENTS:
        
            - text: string containing the FASTA
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
        - flavor: (optional)
            - 0 Protein, L amino acids (default)
            - 1 Protein, D amino acids
            - 2 RNA, no cap
            - 3 RNA, 5' cap
            - 4 RNA, 3' cap
            - 5 RNA, both caps
            - 6 DNA, no cap
            - 7 DNA, 5' cap
            - 8 DNA, 3' cap
            - 9 DNA, both caps
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromFASTA(boost::python::api::object [,bool=True [,int=0]])
    """
def MolFromHELM( text: object, sanitize: bool = True) -> Mol:
    """
    MolFromHELM( text: object, sanitize: bool = True) -> Mol
        Construct a molecule from a HELM string (currently only supports peptides).
        
          ARGUMENTS:
        
            - text: string containing the HELM
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromHELM(boost::python::api::object [,bool=True])
    """
def MolFromMol2Block( molBlock: str, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> Mol:
    """
    MolFromMol2Block( molBlock: str, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> Mol
        Construct a molecule from a Tripos Mol2 block.
        
          NOTE:
            The parser expects the atom-typing scheme used by Corina.
            Atom types from Tripos' dbtranslate are less supported.
            Other atom typing schemes are unlikely to work.
        
          ARGUMENTS:
        
            - mol2Block: string containing the Mol2 block
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - cleanupSubstructures: (optional) toggles standardizing some 
              substructures found in mol2 files.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMol2Block(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
    """
def MolFromMol2File( molFileName: str, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> Mol:
    """
    MolFromMol2File( molFileName: str, sanitize: bool = True, removeHs: bool = True, cleanupSubstructures: bool = True) -> Mol
        Construct a molecule from a Tripos Mol2 file.
        
          NOTE:
            The parser expects the atom-typing scheme used by Corina.
            Atom types from Tripos' dbtranslate are less supported.
            Other atom typing schemes are unlikely to work.
        
          ARGUMENTS:
                                          
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to true.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - cleanupSubstructures: (optional) toggles standardizing some 
              substructures found in mol2 files.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMol2File(char const* [,bool=True [,bool=True [,bool=True]]])
    """
def MolFromMolBlock( molBlock: object, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol:
    """
    MolFromMolBlock( molBlock: object, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol
        Construct a molecule from a Mol block.
        
          ARGUMENTS:
        
            - molBlock: string containing the Mol block
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - strictParsing: (optional) if this is false, the parser is more lax about.
              correctness of the content.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMolBlock(boost::python::api::object [,bool=True [,bool=True [,bool=True]]])

        C++ signature :
            RDKit::ROMol* MolFromMolBlock(boost::python::api::object [,bool=True [,bool=True [,bool=True]]])
    """
def MolFromMolFile( molFileName: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol:
    """
    MolFromMolFile( molFileName: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol
        Construct a molecule from a Mol file.
        
          ARGUMENTS:
        
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to true.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - strictParsing: (optional) if this is false, the parser is more lax about.
              correctness of the content.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMolFile(char const* [,bool=True [,bool=True [,bool=True]]])

        C++ signature :
            RDKit::ROMol* MolFromMolFile(char const* [,bool=True [,bool=True [,bool=True]]])
    """
def MolFromMrvBlock( mrvBlock: object, sanitize: bool = True, removeHs: bool = True) -> Mol:
    """
    MolFromMrvBlock( mrvBlock: object, sanitize: bool = True, removeHs: bool = True) -> Mol
        Construct a molecule from a Marvin (mrv) block.
        
          ARGUMENTS:
        
            - molBlock: string containing the Marvin block
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMrvBlock(boost::python::api::object [,bool=True [,bool=True]])
    """
def MolFromMrvFile( molFileName: str, sanitize: bool = True, removeHs: bool = True) -> Mol:
    """
    MolFromMrvFile( molFileName: str, sanitize: bool = True, removeHs: bool = True) -> Mol
        Construct a molecule from a Marvin (Mrv) file.
        
          ARGUMENTS:
        
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to true.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromMrvFile(char const* [,bool=True [,bool=True]])
    """
def MolFromPDBBlock( molBlock: object, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol:
    """
    MolFromPDBBlock( molBlock: object, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol
        Construct a molecule from a PDB block.
        
          ARGUMENTS:
        
            - molBlock: string containing the PDB block
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - flavor: (optional) 
        
            - proximityBonding: (optional) toggles automatic proximity bonding
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromPDBBlock(boost::python::api::object [,bool=True [,bool=True [,unsigned int=0 [,bool=True]]]])
    """
def MolFromPDBFile( molFileName: str, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol:
    """
    MolFromPDBFile( molFileName: str, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol
        Construct a molecule from a PDB file.
        
          ARGUMENTS:
        
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to true.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
            - flavor: (optional) 
        
            - proximityBonding: (optional) toggles automatic proximity bonding
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromPDBFile(char const* [,bool=True [,bool=True [,unsigned int=0 [,bool=True]]]])
    """
def MolFromPNGFile( filename: str, params: object = None) -> Mol:
    """
    MolFromPNGFile( filename: str, params: object = None) -> Mol
        Construct a molecule from metadata in a PNG file.
        
             ARGUMENTS:
        
               - filename: the PNG filename
        
               - params: used to provide optional parameters for the metadata parsing
        
             RETURNS:
               a Mol object, None on failure.

        C++ signature :
            RDKit::ROMol* MolFromPNGFile(char const* [,boost::python::api::object=None])
    """
def MolFromPNGString( png: object, params: object = None) -> Mol:
    """
    MolFromPNGString( png: object, params: object = None) -> Mol
        Construct a molecule from metadata in a PNG string.
        
             ARGUMENTS:
        
               - png: the PNG string
        
               - params: used to provide optional parameters for the metadata parsing
        
             RETURNS:
               a Mol object, None on failure.
          

        C++ signature :
            RDKit::ROMol* MolFromPNGString(boost::python::api::object [,boost::python::api::object=None])
    """
def MolFromRDKitSVG( molBlock: object, sanitize: bool = True, removeHs: bool = True) -> Mol:
    """
    MolFromRDKitSVG( molBlock: object, sanitize: bool = True, removeHs: bool = True) -> Mol
        Construct a molecule from an RDKit-generate SVG string.
        
          ARGUMENTS:
        
            - svg: string containing the SVG data (must include molecule metadata)
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - removeHs: (optional) toggles removing hydrogens from the molecule.
              This only make sense when sanitization is done.
              Defaults to true.
        
          RETURNS:
        
            a Mol object, None on failure.
        
          NOTE: this functionality should be considered beta.
        
        

        C++ signature :
            RDKit::ROMol* MolFromRDKitSVG(boost::python::api::object [,bool=True [,bool=True]])
    """
def MolFromSequence( text: object, sanitize: bool = True, flavor: int = 0) -> Mol:
    """
    MolFromSequence( text: object, sanitize: bool = True, flavor: int = 0) -> Mol
        Construct a molecule from a sequence string (currently only supports peptides).
        
          ARGUMENTS:
        
            - text: string containing the sequence
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - flavor: (optional)
                - 0 Protein, L amino acids (default)
                - 1 Protein, D amino acids
                - 2 RNA, no cap
                - 3 RNA, 5' cap
                - 4 RNA, 3' cap
                - 5 RNA, both caps
                - 6 DNA, no cap
                - 7 DNA, 5' cap
                - 8 DNA, 3' cap
                - 9 DNA, both caps
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromSequence(boost::python::api::object [,bool=True [,int=0]])
    """
@typing.overload
def MolFromSmarts( SMARTS: object, mergeHs: bool = False, replacements: dict = {}) -> Mol:
    """
    MolFromSmarts( SMARTS: object, mergeHs: bool = False, replacements: dict = {}) -> Mol
        Construct a molecule from a SMARTS string.
        
          ARGUMENTS:
        
            - SMARTS: the smarts string
        
            - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
              atoms.  So, for example, 'C[H]' becomes '[C;!H0]'.
              Defaults to 0.
        
            - replacements: (optional) a dictionary of replacement strings (see below)
              Defaults to {}. See the documentation for MolFromSmiles for an explanation.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromSmarts(boost::python::api::object [,bool=False [,boost::python::dict={}]])

        C++ signature :
            RDKit::ROMol* MolFromSmarts(boost::python::api::object,RDKit::SmartsParserParams)
    """
@typing.overload
def MolFromSmarts( SMARTS: object, params: SmartsParserParams) -> Mol:
    pass
@typing.overload
def MolFromSmiles( SMILES: object, params: SmilesParserParams) -> Mol:
    """
    MolFromSmiles( SMILES: object, params: SmilesParserParams) -> Mol
        Construct a molecule from a SMILES string.
        
             ARGUMENTS:
           
               - SMILES: the smiles string
           
               - params: used to provide optional parameters for the SMILES parsing
           
             RETURNS:
           
               a Mol object, None on failure.
           
        

        C++ signature :
            RDKit::ROMol* MolFromSmiles(boost::python::api::object,RDKit::SmilesParserParams)

        C++ signature :
            RDKit::ROMol* MolFromSmiles(boost::python::api::object [,bool=True [,boost::python::dict={}]])
    """
@typing.overload
def MolFromSmiles( SMILES: object, sanitize: bool = True, replacements: dict = {}) -> Mol:
    pass
def MolFromTPLBlock( tplBlock: object, sanitize: bool = True, skipFirstConf: bool = False) -> Mol:
    """
    MolFromTPLBlock( tplBlock: object, sanitize: bool = True, skipFirstConf: bool = False) -> Mol
        Construct a molecule from a TPL block.
        
          ARGUMENTS:
        
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - skipFirstConf: (optional) skips reading the first conformer.
              Defaults to False.
              This should be set to True when reading TPLs written by 
              the CombiCode.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromTPLBlock(boost::python::api::object [,bool=True [,bool=False]])
    """
def MolFromTPLFile( fileName: str, sanitize: bool = True, skipFirstConf: bool = False) -> Mol:
    """
    MolFromTPLFile( fileName: str, sanitize: bool = True, skipFirstConf: bool = False) -> Mol
        Construct a molecule from a TPL file.
        
          ARGUMENTS:
        
            - fileName: name of the file to read
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
            - skipFirstConf: (optional) skips reading the first conformer.
              Defaults to False.
              This should be set to True when reading TPLs written by 
              the CombiCode.
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromTPLFile(char const* [,bool=True [,bool=False]])
    """
def MolFromXYZBlock( xyzFileName: object) -> Mol:
    """
    MolFromXYZBlock( xyzFileName: object) -> Mol
        Construct a molecule from an XYZ string.
        
          ARGUMENTS:
        
            - xyzBlock: the XYZ data to read
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromXYZBlock(boost::python::api::object)
    """
def MolFromXYZFile( xyzFileName: str) -> Mol:
    """
    MolFromXYZFile( xyzFileName: str) -> Mol
        Construct a molecule from an XYZ file.
        
          ARGUMENTS:
        
            - xyzname: name of the file to read
        
          RETURNS:
        
            a Mol object, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromXYZFile(char const*)
    """
def MolMetadataToPNGFile( mol: Mol, filename: object, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object:
    """
    MolMetadataToPNGFile( mol: Mol, filename: object, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object
        Adds molecular metadata to PNG data read from a file.
        
             ARGUMENTS:
        
               - mol: the molecule
        
               - filename: the PNG filename
        
               - includePkl: include the RDKit's internal binary format in the output
        
               - includeSmiles: include CXSmiles in the output
        
               - includeMol: include CTAB (Mol) in the output
        
             RETURNS:
               the updated PNG data

        C++ signature :
            boost::python::api::object MolMetadataToPNGFile(RDKit::ROMol,boost::python::api::object [,bool=True [,bool=True [,bool=False]]])
    """
def MolMetadataToPNGString( mol: Mol, png: object, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object:
    """
    MolMetadataToPNGString( mol: Mol, png: object, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object
        Adds molecular metadata to a PNG string.
        
             ARGUMENTS:
        
               - mol: the molecule
        
               - png: the PNG string
        
               - includePkl: include the RDKit's internal binary format in the output
        
               - includeSmiles: include CXSmiles in the output
        
               - includeMol: include CTAB (Mol) in the output
        
             RETURNS:
               the updated PNG data

        C++ signature :
            boost::python::api::object MolMetadataToPNGString(RDKit::ROMol,boost::python::api::object [,bool=True [,bool=True [,bool=False]]])
    """
def MolToCMLBlock( mol: Mol, confId: int = -1, kekulize: bool = True) -> str:
    """
    MolToCMLBlock( mol: Mol, confId: int = -1, kekulize: bool = True) -> str
        Writes a CML block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - confId: (optional) selects which conformation to output
            - kekulize: (optional) triggers kekulization of the molecule before it's written
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCMLBlock(RDKit::ROMol [,int=-1 [,bool=True]])
    """
def MolToCMLFile( mol: Mol, filename: str, confId: int = -1, kekulize: bool = True) -> None:
    """
    MolToCMLFile( mol: Mol, filename: str, confId: int = -1, kekulize: bool = True) -> None
        Writes a CML file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: the file to write to
            - confId: (optional) selects which conformation to output
            - kekulize: (optional) triggers kekulization of the molecule before it's written
        
        

        C++ signature :
            void MolToCMLFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1 [,bool=True]])
    """
def MolToCXSmarts( mol: Mol, isomericSmiles: bool = True) -> str:
    """
    MolToCXSmarts( mol: Mol, isomericSmiles: bool = True) -> str
        Returns a SMARTS string for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - isomericSmiles: (optional) include information about stereochemistry in
              the SMARTS.  Defaults to true.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCXSmarts(RDKit::ROMol [,bool=True])
    """
@typing.overload
def MolToCXSmiles( mol: Mol, params: SmilesWriteParams, flags: int = rdmolfiles.CXSmilesFields.CX_ALL) -> str:
    """
    MolToCXSmiles( mol: Mol, params: SmilesWriteParams, flags: int = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL) -> str
        Returns the CXSMILES string for a molecule

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCXSmiles(RDKit::ROMol,RDKit::SmilesWriteParams [,unsigned int=rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL])

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToCXSmiles(RDKit::ROMol [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False [,bool=False]]]]]]])
    """
@typing.overload
def MolToCXSmiles( mol: Mol, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False, doRandom: bool = False) -> str:
    pass
def MolToFASTA( mol: Mol) -> str:
    """
    MolToFASTA( mol: Mol) -> str
        Returns the FASTA string for a molecule
          ARGUMENTS:
        
            - mol: the molecule
        
          NOTE: the molecule should contain monomer information in AtomMonomerInfo structures 
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToFASTA(RDKit::ROMol)
    """
def MolToHELM( mol: Mol) -> str:
    """
    MolToHELM( mol: Mol) -> str
        Returns the HELM string for a molecule
          ARGUMENTS:
        
            - mol: the molecule
        
          NOTE: the molecule should contain monomer information in AtomMonomerInfo structures 
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToHELM(RDKit::ROMol)
    """
def MolToJSON( mol: Mol, params: object = None) -> str:
    """
    MolToJSON( mol: Mol, params: object = None) -> str
        Convert a single molecule to JSON
        
            ARGUMENTS:
              - mol: the molecule to work with
            RETURNS:
              a string
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToJSON(RDKit::ROMol [,boost::python::api::object=None])
    """
def MolToMolBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> str:
    """
    MolToMolBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> str
        Returns a Mol block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written,
              as suggested by the MDL spec.
            - forceV3000 (optional) force generation a V3000 mol block (happens automatically with 
              more than 999 atoms or bonds)
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToMolBlock(RDKit::ROMol [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
def MolToMolFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> None:
    """
    MolToMolFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, forceV3000: bool = False) -> None
        Writes a Mol file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: the file to write to
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written,
              as suggested by the MDL spec.
            - forceV3000 (optional) force generation a V3000 mol block (happens automatically with 
              more than 999 atoms or bonds)
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            void MolToMolFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
def MolToMrvBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> str:
    """
    MolToMrvBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> str
        Returns a Marvin (Mrv) Mol block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written.
            - prettyPrint: (optional) makes the output more human readable.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToMrvBlock(RDKit::ROMol [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
def MolToMrvFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> None:
    """
    MolToMrvFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True, prettyPrint: bool = False) -> None
        Writes a Marvin (MRV) file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: the file to write to
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written.
            - prettyPrint: (optional) makes the output more human readable.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            void MolToMrvFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,int=-1 [,bool=True [,bool=False]]]])
    """
def MolToPDBBlock( mol: Mol, confId: int = -1, flavor: int = 0) -> str:
    """
    MolToPDBBlock( mol: Mol, confId: int = -1, flavor: int = 0) -> str
        Returns a PDB block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - confId: (optional) selects which conformation to output (-1 = default)
            - flavor: (optional) 
                    - flavor & 1 : Write MODEL/ENDMDL lines around each record 
                    - flavor & 2 : Don't write any CONECT records 
                    - flavor & 4 : Write CONECT records in both directions 
                    - flavor & 8 : Don't use multiple CONECTs to encode bond order 
                    - flavor & 16 : Write MASTER record 
                    - flavor & 32 : Write TER record 
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToPDBBlock(RDKit::ROMol [,int=-1 [,unsigned int=0]])
    """
def MolToPDBFile( mol: Mol, filename: str, confId: int = -1, flavor: int = 0) -> None:
    """
    MolToPDBFile( mol: Mol, filename: str, confId: int = -1, flavor: int = 0) -> None
        Writes a PDB file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: name of the file to write
            - confId: (optional) selects which conformation to output (-1 = default)
            - flavor: (optional) 
                    - flavor & 1 : Write MODEL/ENDMDL lines around each record 
                    - flavor & 2 : Don't write any CONECT records 
                    - flavor & 4 : Write CONECT records in both directions 
                    - flavor & 8 : Don't use multiple CONECTs to encode bond order 
                    - flavor & 16 : Write MASTER record 
                    - flavor & 32 : Write TER record 
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            void MolToPDBFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1 [,unsigned int=0]])
    """
def MolToRandomSmilesVect( mol: Mol, numSmiles: int, randomSeed: int = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> list:
    """
    MolToRandomSmilesVect( mol: Mol, numSmiles: int, randomSeed: int = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> list
        returns a list of SMILES generated using the randomSmiles algorithm

        C++ signature :
            boost::python::list MolToRandomSmilesVect(RDKit::ROMol,unsigned int [,unsigned int=0 [,bool=True [,bool=False [,bool=False [,bool=False]]]]])
    """
def MolToSequence( mol: Mol) -> str:
    """
    MolToSequence( mol: Mol) -> str
        Returns the sequence string for a molecule
          ARGUMENTS:
        
            - mol: the molecule
        
          NOTE: the molecule should contain monomer information in AtomMonomerInfo structures 
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSequence(RDKit::ROMol)
    """
def MolToSmarts( mol: Mol, isomericSmiles: bool = True, rootedAtAtom: int = -1) -> str:
    """
    MolToSmarts( mol: Mol, isomericSmiles: bool = True, rootedAtAtom: int = -1) -> str
        Returns a SMARTS string for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - isomericSmiles: (optional) include information about stereochemistry in
              the SMARTS.  Defaults to true.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSmarts(RDKit::ROMol [,bool=True [,int=-1]])
    """
@typing.overload
def MolToSmiles( mol: Mol, params: SmilesWriteParams) -> str:
    """
    MolToSmiles( mol: Mol, params: SmilesWriteParams) -> str
        Returns the canonical SMILES string for a molecule

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSmiles(RDKit::ROMol,RDKit::SmilesWriteParams)

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSmiles(RDKit::ROMol [,bool=True [,bool=False [,int=-1 [,bool=True [,bool=False [,bool=False [,bool=False]]]]]]])
    """
@typing.overload
def MolToSmiles( mol: Mol, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False, doRandom: bool = False) -> str:
    pass
def MolToTPLBlock( mol: Mol, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> str:
    """
    MolToTPLBlock( mol: Mol, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> str
        Returns the Tpl block for a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule
            - partialChargeProp: name of the property to use for partial charges
              Defaults to '_GasteigerCharge'.
            - writeFirstConfTwice: Defaults to False.
              This should be set to True when writing TPLs to be read by 
              the CombiCode.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToTPLBlock(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='_GasteigerCharge' [,bool=False]])
    """
def MolToTPLFile( mol: Mol, fileName: str, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> None:
    """
    MolToTPLFile( mol: Mol, fileName: str, partialChargeProp: str = '_GasteigerCharge', writeFirstConfTwice: bool = False) -> None
        Writes a molecule to a TPL file.
        
          ARGUMENTS:
        
            - mol: the molecule
            - fileName: name of the file to write
            - partialChargeProp: name of the property to use for partial charges
              Defaults to '_GasteigerCharge'.
            - writeFirstConfTwice: Defaults to False.
              This should be set to True when writing TPLs to be read by 
              the CombiCode.
        
        

        C++ signature :
            void MolToTPLFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='_GasteigerCharge' [,bool=False]])
    """
def MolToV3KMolBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> str:
    """
    MolToV3KMolBlock( mol: Mol, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> str
        Returns a V3000 Mol block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written,
              as suggested by the MDL spec.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToV3KMolBlock(RDKit::ROMol [,bool=True [,int=-1 [,bool=True]]])
    """
def MolToV3KMolFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> None:
    """
    MolToV3KMolFile( mol: Mol, filename: str, includeStereo: bool = True, confId: int = -1, kekulize: bool = True) -> None
        Writes a V3000 Mol file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: the file to write to
            - includeStereo: (optional) toggles inclusion of stereochemical
              information in the output
            - confId: (optional) selects which conformation to output (-1 = default)
            - kekulize: (optional) triggers kekulization of the molecule before it's written,
              as suggested by the MDL spec.
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            void MolToV3KMolFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,int=-1 [,bool=True]]])
    """
def MolToXYZBlock( mol: Mol, confId: int = -1) -> str:
    """
    MolToXYZBlock( mol: Mol, confId: int = -1) -> str
        Returns a XYZ block for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - confId: (optional) selects which conformation to output (-1 = default)
        
          RETURNS:
        
            a string
        
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToXYZBlock(RDKit::ROMol [,int=-1])
    """
def MolToXYZFile( mol: Mol, filename: str, confId: int = -1) -> None:
    """
    MolToXYZFile( mol: Mol, filename: str, confId: int = -1) -> None
        Writes a XYZ file for a molecule
          ARGUMENTS:
        
            - mol: the molecule
            - filename: the file to write to
            - confId: (optional) selects which conformation to output (-1 = default)
        
        

        C++ signature :
            void MolToXYZFile(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=-1])
    """
def MolsFromCDXML( cdxml: object, sanitize: bool = True, removeHs: bool = True) -> tuple:
    """
    MolsFromCDXML( cdxml: object, sanitize: bool = True, removeHs: bool = True) -> tuple
        Construct a molecule from a cdxml string.
        
             Note that the CDXML format is large and complex, the RDKit doesn't support
             full functionality, just the base ones required for molecule and
             reaction parsing.
        
             ARGUMENTS:
        
               - filename: the cdxml string
        
               - sanitize: if True, sanitize the molecules [default True]
               - removeHs: if True, convert explicit Hs into implicit Hs. [default True]
        
        
             RETURNS:
               an iterator of parsed Mol objects.

        C++ signature :
            boost::python::tuple MolsFromCDXML(boost::python::api::object [,bool=True [,bool=True]])
    """
def MolsFromCDXMLFile( filename: str, sanitize: bool = True, removeHs: bool = True) -> object:
    """
    MolsFromCDXMLFile( filename: str, sanitize: bool = True, removeHs: bool = True) -> object
        Construct a molecule from a cdxml file.
        
             Note that the CDXML format is large and complex, the RDKit doesn't support
             full functionality, just the base ones required for molecule and
             reaction parsing.
        
             ARGUMENTS:
        
               - filename: the cdxml filename
        
               - sanitize: if True, sanitize the molecules [default True]
               - removeHs: if True, convert explicit Hs into implicit Hs. [default True]
        
             RETURNS:
               an iterator of parsed Mol objects.

        C++ signature :
            boost::python::api::object MolsFromCDXMLFile(char const* [,bool=True [,bool=True]])
    """
def MolsFromPNGFile( filename: str, tag: str = 'rdkitPKL', params: object = None) -> object:
    """
    MolsFromPNGFile( filename: str, tag: str = 'rdkitPKL', params: object = None) -> object
        returns a tuple of molecules constructed from the PNG file

        C++ signature :
            boost::python::api::object MolsFromPNGFile(char const* [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='rdkitPKL' [,boost::python::api::object=None]])
    """
def MolsFromPNGString( png: object, tag: str = 'rdkitPKL', params: object = None) -> tuple:
    """
    MolsFromPNGString( png: object, tag: str = 'rdkitPKL', params: object = None) -> tuple
        returns a tuple of molecules constructed from the PNG string

        C++ signature :
            boost::python::tuple MolsFromPNGString(boost::python::api::object [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='rdkitPKL' [,boost::python::api::object=None]])
    """
def MolsToJSON( mols: object, params: object = None) -> str:
    """
    MolsToJSON( mols: object, params: object = None) -> str
        Convert a set of molecules to JSON
        
            ARGUMENTS:
              - mols: the molecules to work with
            RETURNS:
              a string
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolsToJSON(boost::python::api::object [,boost::python::api::object=None])
    """
def MurckoDecompose( mol: Mol) -> Mol:
    """
    MurckoDecompose( mol: Mol) -> Mol
        Do a Murcko decomposition and return the scaffold

        C++ signature :
            RDKit::ROMol* MurckoDecompose(RDKit::ROMol)
    """
def ParseMolQueryDefFile( fileobj: object, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict:
    """
    ParseMolQueryDefFile( fileobj: object, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict
        reads query definitions from a simply formatted file
        

        C++ signature :
            boost::python::dict ParseMolQueryDefFile(boost::python::api::object {lvalue} [,bool=True [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='\t' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='//' [,unsigned int=0 [,unsigned int=1]]]]])
    """
def PathToSubmol( mol: Mol, path: object, useQuery: bool = False, atomMap: object = None) -> Mol:
    """
    PathToSubmol( mol: Mol, path: object, useQuery: bool = False, atomMap: object = None) -> Mol

        C++ signature :
            RDKit::ROMol* PathToSubmol(RDKit::ROMol,boost::python::api::object {lvalue} [,bool=False [,boost::python::api::object=None]])
    """
@typing.overload
def PatternFingerprint( mol: Mol, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, tautomerFingerprints: bool = False) -> ExplicitBitVect:
    """
    PatternFingerprint( mol: Mol, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, tautomerFingerprints: bool = False) -> ExplicitBitVect
        A fingerprint using SMARTS patterns 
        
          NOTE: This function is experimental. The API or results may change from
            release to release.
        

        C++ signature :
            ExplicitBitVect* PatternFingerprint(RDKit::ROMol [,unsigned int=2048 [,boost::python::list=[] [,ExplicitBitVect*=None [,bool=False]]]])

        C++ signature :
            ExplicitBitVect* PatternFingerprint(RDKit::MolBundle [,unsigned int=2048 [,ExplicitBitVect*=None [,bool=False]]])
    """
@typing.overload
def PatternFingerprint( mol: MolBundle, fpSize: int = 2048, setOnlyBits: ExplicitBitVect = None, tautomerFingerprints: bool = False) -> ExplicitBitVect:
    pass
def RDKFingerprint( mol: Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: object = 0, fromAtoms: object = 0, atomBits: object = None, bitInfo: object = None) -> ExplicitBitVect:
    """
    RDKFingerprint( mol: Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: object = 0, fromAtoms: object = 0, atomBits: object = None, bitInfo: object = None) -> ExplicitBitVect
        Returns an RDKit topological fingerprint for a molecule
        
          Explanation of the algorithm below.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - minPath: (optional) minimum number of bonds to include in the subgraphs
              Defaults to 1.
        
            - maxPath: (optional) maximum number of bonds to include in the subgraphs
              Defaults to 7.
        
            - fpSize: (optional) number of bits in the fingerprint
              Defaults to 2048.
        
            - nBitsPerHash: (optional) number of bits to set per path
              Defaults to 2.
        
            - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
              has explicit Hs.
              Defaults to True.
        
            - tgtDensity: (optional) fold the fingerprint until this minimum density has
              been reached
              Defaults to 0.
        
            - minSize: (optional) the minimum size the fingerprint will be folded to when
              trying to reach tgtDensity
              Defaults to 128.
        
            - branchedPaths: (optional) if set both branched and unbranched paths will be
              used in the fingerprint.
              Defaults to True.
        
            - useBondOrder: (optional) if set both bond orders will be used in the path hashes
              Defaults to True.
        
            - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
              Defaults to empty.
        
            - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
              starting from these atoms will be used.
              Defaults to empty.
        
            - atomBits: (optional) an empty list. If provided, the result will contain a list 
              containing the bits each atom sets.
              Defaults to empty.
        
            - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
              with bits as keys and corresponding bond paths as values.
              Defaults to empty.
        
          RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits
        
          ALGORITHM:
        
           This algorithm functions by find all subgraphs between minPath and maxPath in
           length.  For each subgraph:
        
             1) A hash is calculated.
        
             2) The hash is used to seed a random-number generator
        
             3) _nBitsPerHash_ random numbers are generated and used to set the corresponding
                bits in the fingerprint
        
        
        

        C++ signature :
            ExplicitBitVect* RDKFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=7 [,unsigned int=2048 [,unsigned int=2 [,bool=True [,double=0.0 [,unsigned int=128 [,bool=True [,bool=True [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]]]]]]]])
    """
def ReapplyMolBlockWedging( arg1: Mol) -> None:
    """
    ReapplyMolBlockWedging( arg1: Mol) -> None
        Set the wedging to that which was read from the original
             MolBlock, over-riding anything that was originally there.
        
                  ARGUMENTS:
                
                    - molecule: the molecule to update
                
                
        

        C++ signature :
            void ReapplyMolBlockWedging(RDKit::ROMol {lvalue})
    """
def RemoveAllHs( mol: Mol, sanitize: bool = True) -> Mol:
    """
    RemoveAllHs( mol: Mol, sanitize: bool = True) -> Mol
        Returns a copy of the molecule with all Hs removed.

        C++ signature :
            RDKit::ROMol* RemoveAllHs(RDKit::ROMol [,bool=True])
    """
@typing.overload
def RemoveHs( mol: Mol, implicitOnly: bool = False, updateExplicitCount: bool = False, sanitize: bool = True) -> Mol:
    """
    RemoveHs( mol: Mol, implicitOnly: bool = False, updateExplicitCount: bool = False, sanitize: bool = True) -> Mol
        Removes any hydrogens from the graph of a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - implicitOnly: (optional) if this toggle is set, only implicit Hs will
              be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).
        
            - updateExplicitCount: (optional) if this toggle is set, the explicit H count on atoms with 
              Hs will be updated. Default value is 0 (do not update explicit H count).
        
            - sanitize: (optional) if this toggle is set, the molecule will be sanitized after the Hs
              are removed. Default value is 1 (do sanitize).
        
          RETURNS: a new molecule with the Hs removed
        
          NOTES:
        
            - The original molecule is *not* modified.
            - Hydrogens which aren't connected to a heavy atom will not be
              removed.  This prevents molecules like [H][H] from having
              all atoms removed.
            - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),
              will not be removed.
            - two coordinate Hs, like the central H in C[H-]C, will not be removed
            - Hs connected to dummy atoms will not be removed
            - Hs that are part of the definition of double bond Stereochemistry
              will not be removed
            - Hs that are not connected to anything else will not be removed
        
         

        C++ signature :
            RDKit::ROMol* RemoveHs(RDKit::ROMol [,bool=False [,bool=False [,bool=True]]])

        C++ signature :
            RDKit::ROMol* RemoveHs(RDKit::ROMol,RDKit::MolOps::RemoveHsParameters [,bool=True])
    """
@typing.overload
def RemoveHs( mol: Mol, params: RemoveHsParameters, sanitize: bool = True) -> Mol:
    pass
def RemoveStereochemistry( mol: Mol) -> None:
    """
    RemoveStereochemistry( mol: Mol) -> None
        Removes all stereochemistry info from the molecule.
        
        

        C++ signature :
            void RemoveStereochemistry(RDKit::ROMol {lvalue})
    """
def RenumberAtoms( mol: Mol, newOrder: object) -> Mol:
    """
    RenumberAtoms( mol: Mol, newOrder: object) -> Mol
        Returns a copy of a molecule with renumbered atoms
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - newOrder: the new ordering the atoms (should be numAtoms long)
              for example: if newOrder is [3,2,0,1], then atom 3 in the original 
              molecule will be atom 0 in the new one
        
        
        

        C++ signature :
            RDKit::ROMol* RenumberAtoms(RDKit::ROMol,boost::python::api::object {lvalue})
    """
@typing.overload
def ReplaceCore( mol: Mol, core: Mol, matches: object, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> Mol:
    """
    ReplaceCore( mol: Mol, core: Mol, matches: object, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> Mol
        Removes the core of a molecule and labels the sidechains with dummy atoms based on
        The matches indices given in the matching vector matches.
        Calling:
          ReplaceCore(mol,core,mol.GetSubstructMatch(core))
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - coreQuery: the molecule to be used as a substructure query for recognizing the core
        
            - matches: a matching vector of the type returned by mol.GetSubstructMatch(...)
        
            - replaceDummies: toggles replacement of atoms that match dummies in the query
        
            - labelByIndex: toggles labeling the attachment point dummy atoms with 
              the index of the core atom they're attached to.
        
            - requireDummyMatch: if the molecule has side chains that attach at points not
              flagged with a dummy, it will be rejected (None is returned)
        
          RETURNS: a new molecule with the core removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        EXAMPLES:
        
            >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, ReplaceCore
            >>> mol = MolFromSmiles('C1ONNCC1')
            >>> core = MolFromSmiles('NN')
        
            >>> MolToSmiles(ReplaceCore(mol, core, mol.GetSubstructMatch(core)))
            '[1*]OCCC[2*]'
        
            Since NN is symmetric, we should actually get two matches here if we don't
            uniquify the matches.
        
            >>> [MolToSmiles(ReplaceCore(mol, core, match))
            ...     for match in mol.GetSubstructMatches(core, uniquify=False)]
            ['[1*]OCCC[2*]', '[1*]CCCO[2*]']
        
        

        C++ signature :
            RDKit::ROMol* ReplaceCore(RDKit::ROMol,RDKit::ROMol,boost::python::api::object [,bool=True [,bool=False [,bool=False]]])

        C++ signature :
            RDKit::ROMol* ReplaceCore(RDKit::ROMol,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,bool=False]]]])
    """
@typing.overload
def ReplaceCore( mol: Mol, coreQuery: Mol, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False, useChirality: bool = False) -> Mol:
    pass
def ReplaceSidechains( mol: Mol, coreQuery: Mol, useChirality: bool = False) -> Mol:
    """
    ReplaceSidechains( mol: Mol, coreQuery: Mol, useChirality: bool = False) -> Mol
        Replaces sidechains in a molecule with dummy atoms for their attachment points.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - coreQuery: the molecule to be used as a substructure query for recognizing the core
        
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a new molecule with the sidechains removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'
        
            - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''
        
            - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'
        
        

        C++ signature :
            RDKit::ROMol* ReplaceSidechains(RDKit::ROMol,RDKit::ROMol [,bool=False])
    """
def ReplaceSubstructs( mol: Mol, query: Mol, replacement: Mol, replaceAll: bool = False, replacementConnectionPoint: int = 0, useChirality: bool = False) -> object:
    """
    ReplaceSubstructs( mol: Mol, query: Mol, replacement: Mol, replaceAll: bool = False, replacementConnectionPoint: int = 0, useChirality: bool = False) -> object
        Replaces atoms matching a substructure query in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as a substructure query
        
            - replacement: the molecule to be used as the replacement
        
            - replaceAll: (optional) if this toggle is set, all substructures matching
              the query will be replaced in a single result, otherwise each result will
              contain a separate replacement.
              Default value is False (return multiple replacements)
            - replacementConnectionPoint: (optional) index of the atom in the replacement that
              the bond should be made to.
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a tuple of new molecules with the substructures replaced removed
        
          NOTES:
        
            - The original molecule is *not* modified.
            - A bond is only formed to the remaining atoms, if any, that were bonded 
              to the first atom in the substructure query. (For finer control over
              substructure replacement, consider using ChemicalReaction.)
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - ReplaceSubstructs('CCOC','O[CH3]','NC') -> ('CCNC',)
        
            - ReplaceSubstructs('COCCOC','O[CH3]','NC') -> ('COCCNC','CNCCOC')
        
            - ReplaceSubstructs('COCCOC','O[CH3]','NC',True) -> ('CNCCNC',)
        
            - ReplaceSubstructs('COCCOC','O[CH3]','CN',True,1) -> ('CNCCNC',)
        
            - ReplaceSubstructs('CCOC','[CH3]O','NC') -> ('CC.CN',)
        
        

        C++ signature :
            _object* ReplaceSubstructs(RDKit::ROMol,RDKit::ROMol,RDKit::ROMol [,bool=False [,unsigned int=0 [,bool=False]]])
    """
def SanitizeMol( mol: Mol, sanitizeOps: int = rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors: bool = False) -> SanitizeFlags:
    """
    SanitizeMol( mol: Mol, sanitizeOps: int = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors: bool = False) -> SanitizeFlags
        Kekulize, check valencies, set aromaticity, conjugation and hybridization
        
            - The molecule is modified in place.
        
            - If sanitization fails, an exception will be thrown unless catchErrors is set
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
            - sanitizeOps: (optional) sanitization operations to be carried out
              these should be constructed by or'ing together the
              operations in rdkit.Chem.SanitizeFlags
            - catchErrors: (optional) if provided, instead of raising an exception
              when sanitization fails (the default behavior), the 
              first operation that failed (as defined in rdkit.Chem.SanitizeFlags)
              is returned. Zero is returned on success.
        
        

        C++ signature :
            RDKit::MolOps::SanitizeFlags SanitizeMol(RDKit::ROMol {lvalue} [,unsigned long=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL [,bool=False]])
    """
def SetAllowNontetrahedralChirality( arg1: bool) -> None:
    """
    SetAllowNontetrahedralChirality( arg1: bool) -> None
        toggles recognition of non-tetrahedral chirality from 3D structures

        C++ signature :
            void SetAllowNontetrahedralChirality(bool)
    """
def SetAromaticity( mol: Mol, model: AromaticityModel = rdmolops.AromaticityModel.AROMATICITY_DEFAULT) -> None:
    """
    SetAromaticity( mol: Mol, model: AromaticityModel = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT) -> None
        does aromaticity perception
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - model: the model to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void SetAromaticity(RDKit::ROMol {lvalue} [,RDKit::MolOps::AromaticityModel=rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT])
    """
def SetAtomAlias( atom: Atom, rlabel: str) -> None:
    """
    SetAtomAlias( atom: Atom, rlabel: str) -> None
        Sets the atom's MDL alias text.
        Setting to an empty string clears the alias.

        C++ signature :
            void SetAtomAlias(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def SetAtomRLabel( atom: Atom, rlabel: int) -> None:
    """
    SetAtomRLabel( atom: Atom, rlabel: int) -> None
        Sets the atom's MDL RLabel (this is an integer from 0 to 99).
        Setting to 0 clears the rlabel.

        C++ signature :
            void SetAtomRLabel(RDKit::Atom*,int)
    """
def SetAtomValue( atom: Atom, rlabel: str) -> None:
    """
    SetAtomValue( atom: Atom, rlabel: str) -> None
        Sets the atom's MDL alias text.
        Setting to an empty string clears the alias.

        C++ signature :
            void SetAtomValue(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def SetBondStereoFromDirections( mol: Mol) -> None:
    """
    SetBondStereoFromDirections( mol: Mol) -> None
        Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
        
        

        C++ signature :
            void SetBondStereoFromDirections(RDKit::ROMol {lvalue})
    """
def SetConjugation( mol: Mol) -> None:
    """
    SetConjugation( mol: Mol) -> None
        finds conjugated bonds
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void SetConjugation(RDKit::ROMol {lvalue})
    """
def SetDefaultPickleProperties( arg1: int) -> None:
    """
    SetDefaultPickleProperties( arg1: int) -> None
        Set the current global mol pickler options.

        C++ signature :
            void SetDefaultPickleProperties(unsigned int)
    """
def SetDoubleBondNeighborDirections( mol: Mol, conf: object = None) -> None:
    """
    SetDoubleBondNeighborDirections( mol: Mol, conf: object = None) -> None
        Uses the stereo info on double bonds to set the directions of neighboring single bonds
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
        
        

        C++ signature :
            void SetDoubleBondNeighborDirections(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
def SetGenericQueriesFromProperties( mol: Mol, useAtomLabels: bool = True, useSGroups: bool = True) -> None:
    """
    SetGenericQueriesFromProperties( mol: Mol, useAtomLabels: bool = True, useSGroups: bool = True) -> None
        documentation

        C++ signature :
            void SetGenericQueriesFromProperties(RDKit::ROMol {lvalue} [,bool=True [,bool=True]])
    """
def SetHybridization( mol: Mol) -> None:
    """
    SetHybridization( mol: Mol) -> None
        Assigns hybridization states to atoms
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        

        C++ signature :
            void SetHybridization(RDKit::ROMol {lvalue})
    """
def SetSupplementalSmilesLabel( atom: Atom, label: str) -> None:
    """
    SetSupplementalSmilesLabel( atom: Atom, label: str) -> None
        Sets a supplemental label on an atom that is written to the smiles string.
        
        >>> m = Chem.MolFromSmiles("C")
        >>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')
        >>> Chem.MolToSmiles(m)
        'C<xxx>'
        

        C++ signature :
            void SetSupplementalSmilesLabel(RDKit::Atom*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def SetTerminalAtomCoords( arg1: Mol, arg2: int, arg3: int) -> None:
    """
    SetTerminalAtomCoords( arg1: Mol, arg2: int, arg3: int) -> None
        Sets Cartesian coordinates for a terminal atom.
        
          Useful for growing an atom off a molecule with sensible 
          coordinates based on the geometry of the neighbor.
        
          NOTE: this sets the appropriate coordinates in all of the molecule's conformers 
          ARGUMENTS:
        
            - mol: the molecule the atoms belong to.
            - idx: index of the terminal atom whose coordinates are set.
            - mol: index of the bonded neighbor atom.
        
          RETURNS: Nothing
        
        

        C++ signature :
            void SetTerminalAtomCoords(RDKit::ROMol {lvalue},unsigned int,unsigned int)
    """
def SetUseLegacyStereoPerception( arg1: bool) -> None:
    """
    SetUseLegacyStereoPerception( arg1: bool) -> None
        toggles usage of the legacy stereo perception code

        C++ signature :
            void SetUseLegacyStereoPerception(bool)
    """
def SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier:
    """
    SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier

        C++ signature :
            RDKit::SmilesMolSupplier* SmilesMolSupplierFromText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
    """
def SortMatchesByDegreeOfCoreSubstitution( mol: Mol, core: Mol, matches: object) -> object:
    """
    SortMatchesByDegreeOfCoreSubstitution( mol: Mol, core: Mol, matches: object) -> object
        Postprocesses the results of a mol.GetSubstructMatches(core) call 
        where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
        It returns a copy of matches sorted by decreasing number of non-hydrogen matches 
        to the terminal dummy atoms.
        
          ARGUMENTS:
        
            - mol: the molecule GetSubstructMatches was run on
        
            - core: the molecule used as a substructure query
        
            - matches: the result returned by GetSubstructMatches
        
          RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches 
                   to the terminal dummy atoms
        

        C++ signature :
            _object* SortMatchesByDegreeOfCoreSubstitution(RDKit::ROMol,RDKit::ROMol,boost::python::api::object)
    """
def SplitMolByPDBChainId( mol: Mol, whiteList: object = None, negateList: bool = False) -> dict:
    """
    SplitMolByPDBChainId( mol: Mol, whiteList: object = None, negateList: bool = False) -> dict
        Splits a molecule into pieces based on PDB chain information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by chain id with molecules as the values
        
        

        C++ signature :
            boost::python::dict SplitMolByPDBChainId(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
def SplitMolByPDBResidues( mol: Mol, whiteList: object = None, negateList: bool = False) -> dict:
    """
    SplitMolByPDBResidues( mol: Mol, whiteList: object = None, negateList: bool = False) -> dict
        Splits a molecule into pieces based on PDB residue information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by residue name with molecules as the values
        
        

        C++ signature :
            boost::python::dict SplitMolByPDBResidues(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
def TranslateChiralFlagToStereoGroups( mol: Mol, zeroFlagGroupType: StereoGroupType = rdchem.StereoGroupType.STEREO_AND) -> None:
    """
    TranslateChiralFlagToStereoGroups( mol: Mol, zeroFlagGroupType: StereoGroupType = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND) -> None
        Generate enhanced stereo groups based on the status of the chiral flag property.
        
          Arguments:
           - mol: molecule to be modified
           - zeroFlagGroupType: how to handle non-grouped stereo centers when the
                  chiral flag is set to zero
        
          If the chiral flag is set to a value of 1 then all specified tetrahedral
          chiral centers which are not already in StereoGroups will be added to an
          ABS StereoGroup.
        
          If the chiral flag is set to a value of 0 then all specified tetrahedral
          chiral centers will be added to a StereoGroup of the type zeroFlagGroupType
        
          If there is no chiral flag set (i.e. the property is not present), the
          molecule will not be modified.

        C++ signature :
            void TranslateChiralFlagToStereoGroups(RDKit::ROMol {lvalue} [,RDKit::StereoGroupType=rdkit.Chem.rdchem.StereoGroupType.STEREO_AND])
    """
def UnfoldedRDKFingerprintCountBased( mol: Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: object = 0, fromAtoms: object = 0, atomBits: object = None, bitInfo: object = None) -> ULongSparseIntVect:
    """
    UnfoldedRDKFingerprintCountBased( mol: Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: object = 0, fromAtoms: object = 0, atomBits: object = None, bitInfo: object = None) -> ULongSparseIntVect
        Returns an unfolded count-based version of the RDKit fingerprint for a molecule
        
        ARGUMENTS:
            
                - mol: the molecule to use
            
                - minPath: (optional) minimum number of bonds to include in the subgraphs
                  Defaults to 1.
            
                - maxPath: (optional) maximum number of bonds to include in the subgraphs
                  Defaults to 7.
            
                - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
                  has explicit Hs.
                  Defaults to True.
            
                - branchedPaths: (optional) if set both branched and unbranched paths will be
                  used in the fingerprint.
                  Defaults to True.
            
                - useBondOrder: (optional) if set both bond orders will be used in the path hashes
                  Defaults to True.
            
                - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
                  Defaults to empty.
            
                - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
                  starting from these atoms will be used.
                  Defaults to empty.
            
                - atomBits: (optional) an empty list. If provided, the result will contain a list 
                  containing the bits each atom sets.
                  Defaults to empty.
            
                - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
                  with bits as keys and corresponding bond paths as values.
                  Defaults to empty.
             
             
        

        C++ signature :
            RDKit::SparseIntVect<unsigned long>* UnfoldedRDKFingerprintCountBased(RDKit::ROMol [,unsigned int=1 [,unsigned int=7 [,bool=True [,bool=True [,bool=True [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]]]])
    """
def WedgeBond( arg1: Bond, arg2: int, arg3: Conformer) -> None:
    """
    WedgeBond( arg1: Bond, arg2: int, arg3: Conformer) -> None
        Set the wedging on an individual bond from a molecule.
           The wedging scheme used is that from Mol files.
          ARGUMENTS:
            - bond: the bond to update
            - atom ID: the atom from which to do the wedging
            - conformer: the conformer to use to determine wedge direction
        

        C++ signature :
            void WedgeBond(RDKit::Bond*,unsigned int,RDKit::Conformer const*)
    """
def WedgeMolBonds( mol: Mol, conformer: Conformer, params: BondWedgingParameters = None) -> None:
    """
    WedgeMolBonds( mol: Mol, conformer: Conformer, params: BondWedgingParameters = None) -> None
        Set the wedging on single bonds in a molecule.
           The wedging scheme used is that from Mol files.
        
          ARGUMENTS:
        
            - molecule: the molecule to update
            - conformer: the conformer to use to determine wedge direction
        
        
        

        C++ signature :
            void WedgeMolBonds(RDKit::ROMol {lvalue},RDKit::Conformer const* [,RDKit::Chirality::BondWedgingParameters const*=None])
    """
def _HasSubstructMatchStr( pkl: str, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
    """
    _HasSubstructMatchStr( pkl: str, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool
        This function is included to speed substructure queries from databases, 
        it's probably not of
        general interest.
        
          ARGUMENTS:
            - pkl: a Molecule pickle
        
            - query: a Molecule
        
            - recursionPossible: (optional)
        
            - useChirality: (optional)
        
            - useQueryQueryMatches: use query-query matching logic
        
          RETURNS: True or False
        

        C++ signature :
            bool _HasSubstructMatchStr(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::ROMol [,bool=True [,bool=False [,bool=False]]])
    """
@typing.overload
def molzip( a: Mol, b: Mol, params: MolzipParams = MolzipParams()) -> Mol:
    """
    molzip( a: Mol, b: Mol, params: MolzipParams = MolzipParams()) -> Mol
        zip together two molecules using the given matching parameters

        C++ signature :
            RDKit::ROMol* molzip(RDKit::ROMol,RDKit::ROMol [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7f2c653b0f60>])

        C++ signature :
            RDKit::ROMol* molzip(RDKit::ROMol [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7f2c653ca040>])
    """
@typing.overload
def molzip( a: Mol, params: MolzipParams = MolzipParams()) -> Mol:
    pass
def molzipFragments( mols: object, params: MolzipParams = MolzipParams()) -> Mol:
    """
    molzipFragments( mols: object, params: MolzipParams = MolzipParams()) -> Mol
        zip together multiple molecules from an R group decomposition 
        using the given matching parameters.  The first molecule in the list
        must be the core

        C++ signature :
            RDKit::ROMol* molzipFragments(boost::python::api::object {lvalue} [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7f2c653ca0f0>])
    """
def tossit() -> None:
    """
    tossit() -> None

        C++ signature :
            void tossit()
    """
ADJUST_IGNOREALL = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL
ADJUST_IGNORECHAINS = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
ADJUST_IGNOREDUMMIES = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
ADJUST_IGNOREMAPPED = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED
ADJUST_IGNORENONDUMMIES = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES
ADJUST_IGNORENONE = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE
ADJUST_IGNORERINGS = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS
ALLOW_CHARGE_SEPARATION = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
ALLOW_INCOMPLETE_OCTETS = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
AROMATICITY_CUSTOM = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM
AROMATICITY_DEFAULT = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT
AROMATICITY_MDL = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL
AROMATICITY_RDKIT = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT
AROMATICITY_SIMPLE = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE
AllProps = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
AtomProps = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
BondProps = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
CHI_ALLENE = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
CHI_OCTAHEDRAL = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
CHI_OTHER = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
CHI_SQUAREPLANAR = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
CHI_TETRAHEDRAL = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
CHI_TETRAHEDRAL_CCW = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
CHI_TETRAHEDRAL_CW = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
CHI_TRIGONALBIPYRAMIDAL = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
CHI_UNSPECIFIED = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
COMPOSITE_AND = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
COMPOSITE_OR = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
COMPOSITE_XOR = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
ComputedProps = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
CoordsAsDouble = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
INCHI_AVAILABLE = True
KEKULE_ALL = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
LayeredFingerprint_substructLayers = 7
MolProps = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
NoConformers = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
SANITIZE_ADJUSTHS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
SANITIZE_ALL = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL
SANITIZE_CLEANUP = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
SANITIZE_CLEANUPCHIRALITY = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
SANITIZE_CLEANUP_ORGANOMETALLICS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
SANITIZE_FINDRADICALS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
SANITIZE_KEKULIZE = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
SANITIZE_NONE = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
SANITIZE_PROPERTIES = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
SANITIZE_SETAROMATICITY = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
SANITIZE_SETCONJUGATION = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
SANITIZE_SETHYBRIDIZATION = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
SANITIZE_SYMMRINGS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
STEREO_ABSOLUTE = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
UNCONSTRAINED_ANIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
templDir = '/scratch/toscopa1/src/rdkit/Data/'
