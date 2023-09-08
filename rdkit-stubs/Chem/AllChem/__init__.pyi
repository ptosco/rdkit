""" Import all RDKit chemistry modules

"""
from __future__ import annotations
import rdkit.Chem.AllChem
import typing
from rdkit.Chem.rdFingerprintGenerator import AdditionalOutput
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.Chem.rdmolops import AdjustQueryWhichFlags
from rdkit.Chem.rdmolops import AromaticityModel
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.rdFingerprintGenerator import AtomInvariantsGenerator
from rdkit.Chem.rdchem import AtomKekulizeException
from rdkit.Chem.rdchem import AtomMonomerInfo
from rdkit.Chem.rdchem import AtomMonomerType
from rdkit.Chem.rdchem import AtomPDBResidueInfo
from rdkit.Chem.rdFingerprintGenerator import AtomPairFingerprintOptions
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Chem.rdchem import AtomSanitizeException
from rdkit.Chem.rdchem import AtomValenceException
from rdkit.Chem.rdchem import Bond
from rdkit.Chem.rdchem import BondDir
from rdkit.Chem.rdFingerprintGenerator import BondInvariantsGenerator
from rdkit.Chem.rdchem import BondStereo
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolops import BondWedgingParameters
from rdkit.Chem.rdmolfiles import CXSmilesFields
from rdkit.Chem.rdChemReactions import CartesianProductStrategy
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdchem import CompositeQueryType
from rdkit.Chem.rdchem import Conformer
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdDistGeom import EmbedFailureCauses
from rdkit.Chem.rdDistGeom import EmbedParameters
from rdkit.Chem.rdChemReactions import EnumerateLibrary
from rdkit.Chem.rdChemReactions import EnumerateLibraryBase
from rdkit.Chem.rdChemReactions import EnumerationParams
from rdkit.Chem.rdChemReactions import EnumerationStrategyBase
from rdkit.Chem.rdMolEnumerator import EnumeratorType
from rdkit.Chem.rdChemReactions import EvenSamplePairsStrategy
from rdkit.Chem.rdFingerprintGenerator import FPType
from rdkit.Chem.rdFingerprintGenerator import FingeprintGenerator32
from rdkit.Chem.rdFingerprintGenerator import FingeprintGenerator64
from rdkit.Chem.rdFingerprintGenerator import FingerprintOptions
from rdkit.Chem.rdChemReactions import FingerprintType
from rdkit.Chem.rdchem import FixedMolSizeMolBundle
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdChemicalFeatures import FreeChemicalFeature
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.inchi import InchiReadWriteError
from rdkit.Chem.rdMolInterchange import JSONParseParameters
from rdkit.Chem.rdMolInterchange import JSONWriteParameters
from rdkit.Chem.rdchem import KekulizeException
from rdkit.Chem.rdChemReactions import MOL_SPTR_VECT
from rdkit.Chem.rdmolfiles import MaeMolSupplier
from rdkit.Chem.rdmolfiles import MaeWriter
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdchem import MolBundle
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeatureFactory
from rdkit.Chem.rdMolEnumerator import MolEnumeratorParams
from rdkit.Chem.rdchem import MolSanitizeException
from rdkit.Chem.rdmolops import MolzipLabel
from rdkit.Chem.rdmolops import MolzipParams
from rdkit.Chem.rdFingerprintGenerator import MorganFingerprintOptions
from rdkit.Chem.rdmolfiles import MultithreadedSDMolSupplier
from rdkit.Chem.rdmolfiles import MultithreadedSmilesMolSupplier
from rdkit.Chem.rdMolDescriptors import NumRotatableBondsOptions
from rdkit.Chem.rdMolAlign import O3A
from rdkit.Chem.rdmolfiles import PDBWriter
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem.rdMolDescriptors import PropertyFunctor
from rdkit.Chem.rdchem import PropertyPickleOptions
from rdkit.Chem.rdMolDescriptors import PropertyRangeQuery
from rdkit.Chem.rdMolDescriptors import PythonPropertyFunctor
from rdkit.Chem.rdchem import QueryAtom
from rdkit.Chem.rdchem import QueryBond
from rdkit.Chem.rdFingerprintGenerator import RDKitFingerprintOptions
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdChemReactions import RandomSampleAllBBsStrategy
from rdkit.Chem.rdChemReactions import RandomSampleStrategy
from rdkit.Chem.rdChemReactions import ReactionFingerprintParams
from rdkit.Chem.rdmolops import RemoveHsParameters
from rdkit.Chem.rdchem import ResonanceFlags
from rdkit.Chem.rdchem import ResonanceMolSupplier
from rdkit.Chem.rdchem import ResonanceMolSupplierCallback
from rdkit.Chem.rdchem import RingInfo
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem.rdChemReactions import SanitizeFlags
from rdkit.Chem.rdmolfiles import SmartsParserParams
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from rdkit.Chem.rdmolfiles import SmilesParserParams
from rdkit.Chem.rdmolfiles import SmilesWriteParams
from rdkit.Chem.rdmolfiles import SmilesWriter
from rdkit.Chem.rdmolops import StereoBondThresholds
from rdkit.Chem.rdchem import StereoDescriptor
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions
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
from rdkit.Chem.rdFingerprintGenerator import TopologicalTorsionFingerprintOptions
from rdkit.Chem.rdDepictor import UsingCoordGen
from rdkit.Chem.rdChemReactions import VectMolVect
from rdkit.Chem.rdChemReactions import VectSizeT
from rdkit.Chem.rdChemReactions import VectorOfStringVectors
import numpy
import rdkit.Chem.inchi
import rdkit.Chem.rdCIPLabeler
import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdCoordGen
import rdkit.Chem.rdDistGeom
import rdkit.Chem.rdFingerprintGenerator
import rdkit.Chem.rdMolInterchange
import rdkit.Chem.rdchem
import rdkit.Chem.rdinchi
import rdkit.Chem.rdmolfiles
import rdkit.Chem.rdmolops
import rdkit.DataStructs
import rdkit.ForceField
import rdkit.Geometry.rdGeometry
import rdkit.RDConfig
import rdkit.RDLogger
import rdkit.rdBase
import sys
import warnings
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AAtomQueryAtom",
    "ADJUST_IGNOREALL",
    "ADJUST_IGNORECHAINS",
    "ADJUST_IGNOREDUMMIES",
    "ADJUST_IGNOREMAPPED",
    "ADJUST_IGNORENONDUMMIES",
    "ADJUST_IGNORENONE",
    "ADJUST_IGNORERINGS",
    "AHAtomQueryAtom",
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
    "AddRingSystemTemplates",
    "AddWavyBondsForStereoAny",
    "AdditionalOutput",
    "AdjustQueryParameters",
    "AdjustQueryProperties",
    "AdjustQueryPropertiesWithGenericGroups",
    "AdjustQueryWhichFlags",
    "AlignMol",
    "AlignMolConformers",
    "AllProps",
    "AromaticityModel",
    "AssignAtomChiralTagsFromMolParity",
    "AssignAtomChiralTagsFromStructure",
    "AssignBondOrdersFromTemplate",
    "AssignCIPLabels",
    "AssignChiralTypesFromBondDirs",
    "AssignRadicals",
    "AssignStereochemistry",
    "AssignStereochemistryFrom3D",
    "Atom",
    "AtomFromSmarts",
    "AtomFromSmiles",
    "AtomInvariantsGenerator",
    "AtomKekulizeException",
    "AtomMonomerInfo",
    "AtomMonomerType",
    "AtomNumEqualsQueryAtom",
    "AtomNumGreaterQueryAtom",
    "AtomNumLessQueryAtom",
    "AtomPDBResidueInfo",
    "AtomPairFP",
    "AtomPairFingerprintOptions",
    "AtomPairsParameters",
    "AtomProps",
    "AtomSanitizeException",
    "AtomValenceException",
    "BAD_DOUBLE_BOND_STEREO",
    "BCUT2D",
    "Bond",
    "BondDir",
    "BondFromSmarts",
    "BondFromSmiles",
    "BondInvariantsGenerator",
    "BondProps",
    "BondStereo",
    "BondType",
    "BondWedgingParameters",
    "BuildFeatureFactory",
    "BuildFeatureFactoryFromString",
    "CHECK_CHIRAL_CENTERS",
    "CHECK_TETRAHEDRAL_CENTERS",
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
    "CalcAUTOCORR2D",
    "CalcAUTOCORR3D",
    "CalcAsphericity",
    "CalcChi0n",
    "CalcChi0v",
    "CalcChi1n",
    "CalcChi1v",
    "CalcChi2n",
    "CalcChi2v",
    "CalcChi3n",
    "CalcChi3v",
    "CalcChi4n",
    "CalcChi4v",
    "CalcChiNn",
    "CalcChiNv",
    "CalcCoulombMat",
    "CalcCrippenDescriptors",
    "CalcEEMcharges",
    "CalcEccentricity",
    "CalcExactMolWt",
    "CalcFractionCSP3",
    "CalcGETAWAY",
    "CalcHallKierAlpha",
    "CalcInertialShapeFactor",
    "CalcKappa1",
    "CalcKappa2",
    "CalcKappa3",
    "CalcLabuteASA",
    "CalcMORSE",
    "CalcMolFormula",
    "CalcNPR1",
    "CalcNPR2",
    "CalcNumAliphaticCarbocycles",
    "CalcNumAliphaticHeterocycles",
    "CalcNumAliphaticRings",
    "CalcNumAmideBonds",
    "CalcNumAromaticCarbocycles",
    "CalcNumAromaticHeterocycles",
    "CalcNumAromaticRings",
    "CalcNumAtomStereoCenters",
    "CalcNumAtoms",
    "CalcNumBridgeheadAtoms",
    "CalcNumHBA",
    "CalcNumHBD",
    "CalcNumHeavyAtoms",
    "CalcNumHeteroatoms",
    "CalcNumHeterocycles",
    "CalcNumLipinskiHBA",
    "CalcNumLipinskiHBD",
    "CalcNumRings",
    "CalcNumRotatableBonds",
    "CalcNumSaturatedCarbocycles",
    "CalcNumSaturatedHeterocycles",
    "CalcNumSaturatedRings",
    "CalcNumSpiroAtoms",
    "CalcNumUnspecifiedAtomStereoCenters",
    "CalcOxidationNumbers",
    "CalcPBF",
    "CalcPMI1",
    "CalcPMI2",
    "CalcPMI3",
    "CalcPhi",
    "CalcRDF",
    "CalcRMS",
    "CalcRadiusOfGyration",
    "CalcSpherocityIndex",
    "CalcTPSA",
    "CalcWHIM",
    "CanonSmiles",
    "CanonicalRankAtoms",
    "CanonicalRankAtomsInFragment",
    "CanonicalizeConformer",
    "CanonicalizeEnhancedStereo",
    "CanonicalizeMol",
    "CartesianProductStrategy",
    "ChemicalReaction",
    "ChiralType",
    "Cleanup",
    "CleanupOrganometallics",
    "ClearMolSubstanceGroups",
    "CombineMols",
    "CompositeQueryType",
    "Compute2DCoords",
    "Compute2DCoordsForReaction",
    "Compute2DCoordsMimicDistmat",
    "ComputeCanonicalTransform",
    "ComputeCentroid",
    "ComputeConfBox",
    "ComputeConfDimsAndOffset",
    "ComputeGasteigerCharges",
    "ComputeMolShape",
    "ComputeMolVolume",
    "ComputePrincipalAxesAndMoments",
    "ComputePrincipalAxesAndMomentsFromGyrationMatrix",
    "ComputeUnionBox",
    "ComputedProps",
    "Conformer",
    "ConstrainedEmbed",
    "ConvertGenericQueriesToSubstanceGroups",
    "CoordsAsDouble",
    "CreateAtomBoolPropertyList",
    "CreateAtomDoublePropertyList",
    "CreateAtomIntPropertyList",
    "CreateAtomStringPropertyList",
    "CreateDifferenceFingerprintForReaction",
    "CreateMolDataSubstanceGroup",
    "CreateMolSubstanceGroup",
    "CreateStereoGroup",
    "CreateStructuralFingerprintForReaction",
    "CustomProp_VSA_",
    "DataStructs",
    "DativeBondsToHaptic",
    "DeleteSubstructs",
    "DetectBondStereoChemistry",
    "DetectBondStereochemistry",
    "DetectChemistryProblems",
    "ETDG",
    "ETKDG",
    "ETKDGv2",
    "ETKDGv3",
    "ETK_MINIMIZATION",
    "EditableMol",
    "EmbedFailureCauses",
    "EmbedMolecule",
    "EmbedMultipleConfs",
    "EmbedParameters",
    "EncodeShape",
    "Enumerate",
    "EnumerateLibrary",
    "EnumerateLibraryBase",
    "EnumerateLibraryCanSerialize",
    "EnumerateLibraryFromReaction",
    "EnumerateStereoisomers",
    "EnumerationParams",
    "EnumerationStrategyBase",
    "EnumeratorType",
    "EvenSamplePairsStrategy",
    "ExplicitDegreeEqualsQueryAtom",
    "ExplicitDegreeGreaterQueryAtom",
    "ExplicitDegreeLessQueryAtom",
    "ExplicitValenceEqualsQueryAtom",
    "ExplicitValenceGreaterQueryAtom",
    "ExplicitValenceLessQueryAtom",
    "FINAL_CENTER_IN_VOLUME",
    "FINAL_CHIRAL_BOUNDS",
    "FIRST_MINIMIZATION",
    "FPType",
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
    "FingeprintGenerator32",
    "FingeprintGenerator64",
    "FingerprintOptions",
    "FingerprintType",
    "FixedMolSizeMolBundle",
    "ForceField",
    "FormalChargeEqualsQueryAtom",
    "FormalChargeGreaterQueryAtom",
    "FormalChargeLessQueryAtom",
    "ForwardSDMolSupplier",
    "ForwardStereoGroupIds",
    "FragmentOnBRICSBonds",
    "FragmentOnBonds",
    "FragmentOnSomeBonds",
    "FreeChemicalFeature",
    "GenerateDepictionMatching2DStructure",
    "GenerateDepictionMatching3DStructure",
    "GenerateErGFingerprintForReducedGraph",
    "GenerateMolExtendedReducedGraph",
    "Get3DDistanceMatrix",
    "GetAdjacencyMatrix",
    "GetAlignmentTransform",
    "GetAllowNontetrahedralChirality",
    "GetAngleDeg",
    "GetAngleRad",
    "GetAtomAlias",
    "GetAtomFeatures",
    "GetAtomMatch",
    "GetAtomPairAtomCode",
    "GetAtomPairAtomInvGen",
    "GetAtomPairCode",
    "GetAtomPairFingerprint",
    "GetAtomPairGenerator",
    "GetAtomRLabel",
    "GetAtomValue",
    "GetBestAlignmentTransform",
    "GetBestRMS",
    "GetBondLength",
    "GetChemDrawRxnAdjustParams",
    "GetConformerRMS",
    "GetConformerRMSMatrix",
    "GetConnectivityInvariants",
    "GetCountFPs",
    "GetCrippenO3A",
    "GetCrippenO3AForProbeConfs",
    "GetDefaultAdjustParams",
    "GetDefaultPickleProperties",
    "GetDihedralDeg",
    "GetDihedralRad",
    "GetDistanceMatrix",
    "GetErGFingerprint",
    "GetExperimentalTorsions",
    "GetFPs",
    "GetFeatureInvariants",
    "GetFormalCharge",
    "GetHashedAtomPairFingerprint",
    "GetHashedAtomPairFingerprintAsBitVect",
    "GetHashedMorganFingerprint",
    "GetHashedTopologicalTorsionFingerprint",
    "GetHashedTopologicalTorsionFingerprintAsBitVect",
    "GetMACCSKeysFingerprint",
    "GetMolFrags",
    "GetMolSubstanceGroupWithIdx",
    "GetMolSubstanceGroups",
    "GetMoleculeBoundsMatrix",
    "GetMorganAtomInvGen",
    "GetMorganBondInvGen",
    "GetMorganFeatureAtomInvGen",
    "GetMorganFingerprint",
    "GetMorganFingerprintAsBitVect",
    "GetMorganGenerator",
    "GetMostSubstitutedCoreMatch",
    "GetO3A",
    "GetO3AForProbeConfs",
    "GetPeriodicTable",
    "GetPreferCoordGen",
    "GetRDKitAtomInvGen",
    "GetRDKitFPGenerator",
    "GetSSSR",
    "GetShortestPath",
    "GetSparseCountFPs",
    "GetSparseFPs",
    "GetSupplementalSmilesLabel",
    "GetSymmSSSR",
    "GetTopologicalTorsionFingerprint",
    "GetTopologicalTorsionGenerator",
    "GetUFFAngleBendParams",
    "GetUFFBondStretchParams",
    "GetUFFInversionParams",
    "GetUFFTorsionParams",
    "GetUFFVdWParams",
    "GetUSR",
    "GetUSRCAT",
    "GetUSRDistributions",
    "GetUSRDistributionsFromPoints",
    "GetUSRFromDistributions",
    "GetUSRScore",
    "GetUseLegacyStereoPerception",
    "HCountEqualsQueryAtom",
    "HCountGreaterQueryAtom",
    "HCountLessQueryAtom",
    "HapticBondsToDative",
    "HasAgentTemplateSubstructMatch",
    "HasBitVectPropWithValueQueryAtom",
    "HasBoolPropWithValueQueryAtom",
    "HasBoolPropWithValueQueryBond",
    "HasChiralTagQueryAtom",
    "HasDoublePropWithValueQueryAtom",
    "HasDoublePropWithValueQueryBond",
    "HasIntPropWithValueQueryAtom",
    "HasIntPropWithValueQueryBond",
    "HasProductTemplateSubstructMatch",
    "HasPropQueryAtom",
    "HasPropQueryBond",
    "HasReactantTemplateSubstructMatch",
    "HasReactionAtomMapping",
    "HasReactionSubstructMatch",
    "HasStringPropWithValueQueryAtom",
    "HasStringPropWithValueQueryBond",
    "HybridizationEqualsQueryAtom",
    "HybridizationGreaterQueryAtom",
    "HybridizationLessQueryAtom",
    "HybridizationType",
    "INCHI_AVAILABLE",
    "INITIAL_COORDS",
    "InNRingsEqualsQueryAtom",
    "InNRingsGreaterQueryAtom",
    "InNRingsLessQueryAtom",
    "InchiReadWriteError",
    "InchiToInchiKey",
    "IsAliphaticQueryAtom",
    "IsAromaticQueryAtom",
    "IsBridgeheadQueryAtom",
    "IsCoordGenSupportAvailable",
    "IsInRingQueryAtom",
    "IsReactionTemplateMoleculeAgent",
    "IsUnsaturatedQueryAtom",
    "IsotopeEqualsQueryAtom",
    "IsotopeGreaterQueryAtom",
    "IsotopeLessQueryAtom",
    "JSONParseParameters",
    "JSONToMols",
    "JSONWriteParameters",
    "KDG",
    "KEKULE_ALL",
    "Kekulize",
    "KekulizeException",
    "KekulizeIfPossible",
    "LINEAR_DOUBLE_BOND",
    "LayeredFingerprint",
    "LayeredFingerprint_substructLayers",
    "LoadDefaultRingSystemTemplates",
    "MAtomQueryAtom",
    "MCFF_GetFeaturesForMol",
    "MHAtomQueryAtom",
    "MINIMIZE_FOURTH_DIMENSION",
    "MMFFGetMoleculeForceField",
    "MMFFGetMoleculeProperties",
    "MMFFHasAllMoleculeParams",
    "MMFFOptimizeMolecule",
    "MMFFOptimizeMoleculeConfs",
    "MMFFSanitizeMolecule",
    "MOL_SPTR_VECT",
    "MQNs_",
    "MaeMolSupplier",
    "MaeWriter",
    "MakePropertyRangeQuery",
    "MassEqualsQueryAtom",
    "MassGreaterQueryAtom",
    "MassLessQueryAtom",
    "MatchOnlyAtRgroupsAdjustParams",
    "MergeQueryHs",
    "MetadataFromPNGFile",
    "MetadataFromPNGString",
    "MinRingSizeEqualsQueryAtom",
    "MinRingSizeGreaterQueryAtom",
    "MinRingSizeLessQueryAtom",
    "MissingChiralTagQueryAtom",
    "Mol",
    "MolAddRecursiveQueries",
    "MolBlockToInchi",
    "MolBlockToInchiAndAuxInfo",
    "MolBundle",
    "MolBundleCanSerialize",
    "MolChemicalFeature",
    "MolChemicalFeatureFactory",
    "MolEnumeratorParams",
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
    "MolFromQuerySLN",
    "MolFromRDKitSVG",
    "MolFromSLN",
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
    "MorganFP",
    "MorganFingerprintOptions",
    "MrvBlockIsReaction",
    "MrvFileIsReaction",
    "MultithreadedSDMolSupplier",
    "MultithreadedSmilesMolSupplier",
    "MurckoDecompose",
    "NoConformers",
    "NoProps",
    "NonHydrogenDegreeEqualsQueryAtom",
    "NonHydrogenDegreeGreaterQueryAtom",
    "NonHydrogenDegreeLessQueryAtom",
    "NormalizeDepiction",
    "NumAliphaticHeteroatomNeighborsEqualsQueryAtom",
    "NumAliphaticHeteroatomNeighborsGreaterQueryAtom",
    "NumAliphaticHeteroatomNeighborsLessQueryAtom",
    "NumHeteroatomNeighborsEqualsQueryAtom",
    "NumHeteroatomNeighborsGreaterQueryAtom",
    "NumHeteroatomNeighborsLessQueryAtom",
    "NumRadicalElectronsEqualsQueryAtom",
    "NumRadicalElectronsGreaterQueryAtom",
    "NumRadicalElectronsLessQueryAtom",
    "NumRotatableBondsOptions",
    "O3A",
    "OptimizeMolecule",
    "OptimizeMoleculeConfs",
    "PDBWriter",
    "PEOE_VSA_",
    "ParseMolQueryDefFile",
    "PathToSubmol",
    "PatternFingerprint",
    "PeriodicTable",
    "PreprocessReaction",
    "PrivateProps",
    "Properties",
    "PropertyFunctor",
    "PropertyPickleOptions",
    "PropertyRangeQuery",
    "PythonPropertyFunctor",
    "QAtomQueryAtom",
    "QHAtomQueryAtom",
    "QueryAtom",
    "QueryAtomData",
    "QueryBond",
    "QuickSmartsMatch",
    "RDConfig",
    "RDKFingerprint",
    "RDKitFP",
    "RDKitFingerprintOptions",
    "RWMol",
    "RandomSampleAllBBsStrategy",
    "RandomSampleStrategy",
    "RandomTransform",
    "ReactionFingerprintParams",
    "ReactionFromMolecule",
    "ReactionFromMrvBlock",
    "ReactionFromMrvFile",
    "ReactionFromPNGFile",
    "ReactionFromPNGString",
    "ReactionFromRxnBlock",
    "ReactionFromRxnFile",
    "ReactionFromSmarts",
    "ReactionMetadataToPNGFile",
    "ReactionMetadataToPNGString",
    "ReactionToMolecule",
    "ReactionToMrvBlock",
    "ReactionToMrvFile",
    "ReactionToRxnBlock",
    "ReactionToSmarts",
    "ReactionToSmiles",
    "ReactionToV3KRxnBlock",
    "ReactionsFromCDXMLBlock",
    "ReactionsFromCDXMLFile",
    "ReapplyMolBlockWedging",
    "ReduceProductToSideChains",
    "RemoveAllHs",
    "RemoveHs",
    "RemoveHsParameters",
    "RemoveMappingNumbersFromReactions",
    "RemoveStereochemistry",
    "RenumberAtoms",
    "ReplaceCore",
    "ReplaceSidechains",
    "ReplaceSubstructs",
    "ResonanceFlags",
    "ResonanceMolSupplier",
    "ResonanceMolSupplierCallback",
    "RingBondCountEqualsQueryAtom",
    "RingBondCountGreaterQueryAtom",
    "RingBondCountLessQueryAtom",
    "RingInfo",
    "SANITIZE_ADJUSTHS",
    "SANITIZE_ADJUST_REACTANTS",
    "SANITIZE_ALL",
    "SANITIZE_ATOM_MAPS",
    "SANITIZE_CLEANUP",
    "SANITIZE_CLEANUPCHIRALITY",
    "SANITIZE_CLEANUP_ORGANOMETALLICS",
    "SANITIZE_FINDRADICALS",
    "SANITIZE_KEKULIZE",
    "SANITIZE_MERGEHS",
    "SANITIZE_NONE",
    "SANITIZE_PROPERTIES",
    "SANITIZE_RGROUP_NAMES",
    "SANITIZE_SETAROMATICITY",
    "SANITIZE_SETCONJUGATION",
    "SANITIZE_SETHYBRIDIZATION",
    "SANITIZE_SYMMRINGS",
    "SDMolSupplier",
    "SDWriter",
    "SMR_VSA_",
    "STEREO_ABSOLUTE",
    "STEREO_AND",
    "STEREO_OR",
    "SanitizeFlags",
    "SanitizeMol",
    "SanitizeRxn",
    "SetAllowNontetrahedralChirality",
    "SetAngleDeg",
    "SetAngleRad",
    "SetAromaticity",
    "SetAtomAlias",
    "SetAtomRLabel",
    "SetAtomValue",
    "SetBondLength",
    "SetBondStereoFromDirections",
    "SetConjugation",
    "SetDefaultPickleProperties",
    "SetDihedralDeg",
    "SetDihedralRad",
    "SetDoubleBondNeighborDirections",
    "SetGenericQueriesFromProperties",
    "SetHybridization",
    "SetPreferCoordGen",
    "SetRingSystemTemplates",
    "SetSupplementalSmilesLabel",
    "SetTerminalAtomCoords",
    "SetUseLegacyStereoPerception",
    "ShapeProtrudeDist",
    "ShapeTanimotoDist",
    "ShapeTverskyIndex",
    "SlogP_VSA_",
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
    "StereoEnumerationOptions",
    "StereoGroup",
    "StereoGroupType",
    "StereoGroup_vect",
    "StereoInfo",
    "StereoSpecified",
    "StereoType",
    "StraightenDepiction",
    "SubstanceGroup",
    "SubstanceGroupAttach",
    "SubstanceGroupCState",
    "SubstanceGroup_VECT",
    "SubstructMatchParameters",
    "SupplierFromFilename",
    "TDTMolSupplier",
    "TDTWriter",
    "TopologicalTorsionFP",
    "TopologicalTorsionFingerprintOptions",
    "TotalDegreeEqualsQueryAtom",
    "TotalDegreeGreaterQueryAtom",
    "TotalDegreeLessQueryAtom",
    "TotalValenceEqualsQueryAtom",
    "TotalValenceGreaterQueryAtom",
    "TotalValenceLessQueryAtom",
    "TransformConformer",
    "TransformMol",
    "TranslateChiralFlagToStereoGroups",
    "UFFGetMoleculeForceField",
    "UFFHasAllMoleculeParams",
    "UFFOptimizeMolecule",
    "UFFOptimizeMoleculeConfs",
    "UNCONSTRAINED_ANIONS",
    "UNCONSTRAINED_CATIONS",
    "UnfoldedRDKFingerprintCountBased",
    "UpdateProductsStereochemistry",
    "UsingCoordGen",
    "VectMolVect",
    "VectSizeT",
    "VectorOfStringVectors",
    "WedgeBond",
    "WedgeMolBonds",
    "XAtomQueryAtom",
    "XHAtomQueryAtom",
    "inchi",
    "logger",
    "molzip",
    "molzipFragments",
    "namedtuple",
    "numpy",
    "rdBase",
    "rdCIPLabeler",
    "rdCoordGen",
    "rdGeometry",
    "rdMolInterchange",
    "rdchem",
    "rdinchi",
    "rdmolfiles",
    "rdmolops",
    "srETKDGv3",
    "sys",
    "templDir",
    "tossit",
    "warnings"
]


def AAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    AAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when AAtom is True.

        C++ signature :
            RDKit::QueryAtom* AAtomQueryAtom([ bool=False])
    """
def AHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    AHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when AHAtom is True.

        C++ signature :
            RDKit::QueryAtom* AHAtomQueryAtom([ bool=False])
    """
def AddHs( mol: Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: AtomPairsParameters = None, addResidueInfo: bool = False) -> Mol:
    """
    AddHs( mol: Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: AtomPairsParameters = None, addResidueInfo: bool = False) -> Mol
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
def AddMetadataToPNGFile( metadata: dict, filename: AtomPairsParameters) -> object:
    """
    AddMetadataToPNGFile( metadata: dict, filename: AtomPairsParameters) -> object
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
def AddMetadataToPNGString( metadata: dict, png: AtomPairsParameters) -> object:
    """
    AddMetadataToPNGString( metadata: dict, png: AtomPairsParameters) -> object
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
def AddRingSystemTemplates( templatePath: str) -> None:
    """
    AddRingSystemTemplates( templatePath: str) -> None
        Adds the ring system templates from the specified file to be used in 2D coordinate generation. If there are duplicates, the most recently added template will be used. Each template must be a single line in the file represented using CXSMILES, and the structure should be a single ring system. Throws a DepictException if any templates are invalid.

        C++ signature :
            void AddRingSystemTemplates(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
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
def AdjustQueryProperties( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    AdjustQueryProperties( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Returns a new molecule where the query properties of atoms have been modified.

        C++ signature :
            RDKit::ROMol* AdjustQueryProperties(RDKit::ROMol [,boost::python::api::object=None])
    """
def AdjustQueryPropertiesWithGenericGroups( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    AdjustQueryPropertiesWithGenericGroups( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Returns a new molecule where the query properties of atoms have been modified and generic group queries have been prepared.

        C++ signature :
            RDKit::ROMol* AdjustQueryPropertiesWithGenericGroups(RDKit::ROMol [,boost::python::api::object=None])
    """
def AlignMol( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, atomMap: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> float:
    """
    AlignMol( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, atomMap: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> float
        Optimally (minimum RMSD) align a molecule to another molecule
             
              The 3D transformation required to align the specied conformation in the probe molecule
              to a specified conformation in the reference molecule is computed so that the root mean
              squared distance between a specified set of atoms is minimized. 
              This transform is then applied to the specified conformation in the probe molecule
             
             ARGUMENTS
              - prbMol    molecule that is to be aligned
              - refMol    molecule used as the reference for the alignment
              - prbCid    ID of the conformation in the probe to be used 
                               for the alignment (defaults to first conformation)
              - refCid    ID of the conformation in the ref molecule to which 
                               the alignment is computed (defaults to first conformation)
              - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                               used to compute the alignments. If this mapping is 
                               not specified an attempt is made to generate one by
                               substructure matching
              - weights   Optionally specify weights for each of the atom pairs
              - reflect   if true reflect the conformation of the probe molecule
              - maxIters  maximum number of iterations used in minimizing the RMSD
               
              RETURNS
              RMSD value
            
        

        C++ signature :
            double AlignMol(RDKit::ROMol {lvalue},RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]])
    """
def AlignMolConformers( mol: Mol, atomIds: AtomPairsParameters = [], confIds: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50, RMSlist: AtomPairsParameters = None) -> None:
    """
    AlignMolConformers( mol: Mol, atomIds: AtomPairsParameters = [], confIds: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50, RMSlist: AtomPairsParameters = None) -> None
        Align conformations in a molecule to each other
             
              The first conformation in the molecule is used as the reference
             
             ARGUMENTS
              - mol          molecule of interest
              - atomIds      List of atom ids to use a points for alingment - defaults to all atoms
              - confIds      Ids of conformations to align - defaults to all conformers 
              - weights      Optionally specify weights for each of the atom pairs
              - reflect      if true reflect the conformation of the probe molecule
              - maxIters     maximum number of iterations used in minimizing the RMSD
              - RMSlist      if provided, fills in the RMS values between the reference
                         conformation and the other aligned conformations
               
            
        

        C++ signature :
            void AlignMolConformers(RDKit::ROMol {lvalue} [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50 [,boost::python::api::object=None]]]]]])
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
def AssignCIPLabels( mol: Mol, atomsToLabel: AtomPairsParameters = None, bondsToLabel: AtomPairsParameters = None, maxRecursiveIterations: int = 0) -> None:
    """
    AssignCIPLabels( mol: Mol, atomsToLabel: AtomPairsParameters = None, bondsToLabel: AtomPairsParameters = None, maxRecursiveIterations: int = 0) -> None
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
def AtomNumEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* AtomNumEqualsQueryAtom(int [,bool=False])
    """
def AtomNumGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* AtomNumGreaterQueryAtom(int [,bool=False])
    """
def AtomNumLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* AtomNumLessQueryAtom(int [,bool=False])
    """
@typing.overload
def BCUT2D( mol: Mol) -> list:
    """
    BCUT2D( mol: Mol) -> list
        Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999Diagonal elements are (currently) atomic mass, gasteiger charge,crippen logP and crippen MRReturns the 2D BCUT2D descriptors vector as described in
        returns [mass eigen value high, mass eigen value low,
                 gasteiger charge eigenvalue high, gasteiger charge low,
                 crippen lowgp  eigenvalue high, crippen lowgp  low,
                 crippen mr eigenvalue high, crippen mr low]
        

        C++ signature :
            boost::python::list BCUT2D(RDKit::ROMol)

        C++ signature :
            std::pair<double, double> BCUT2D(RDKit::ROMol,boost::python::list)

        C++ signature :
            std::pair<double, double> BCUT2D(RDKit::ROMol,boost::python::tuple)

        C++ signature :
            std::pair<double, double> BCUT2D(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
@typing.overload
def BCUT2D( mol: Mol, atom_props: list) -> tuple:
    pass
@typing.overload
def BCUT2D( mol: Mol, atom_props: tuple) -> tuple:
    pass
@typing.overload
def BCUT2D( mol: Mol, atom_propname: str) -> tuple:
    pass
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
def CalcAUTOCORR2D( mol: Mol, CustomAtomProperty: str = '') -> list:
    """
    CalcAUTOCORR2D( mol: Mol, CustomAtomProperty: str = '') -> list
        Returns 2D Autocorrelation descriptors vector

        C++ signature :
            boost::python::list CalcAUTOCORR2D(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
    """
def CalcAUTOCORR3D( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
    CalcAUTOCORR3D( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list
        Returns 3D Autocorrelation descriptors vector

        C++ signature :
            boost::python::list CalcAUTOCORR3D(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']])
    """
def CalcAsphericity( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcAsphericity( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcAsphericity(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcChi0n( mol: Mol, force: bool = False) -> float:
    """
    CalcChi0n( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi0n(RDKit::ROMol [,bool=False])
    """
def CalcChi0v( mol: Mol, force: bool = False) -> float:
    """
    CalcChi0v( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi0v(RDKit::ROMol [,bool=False])
    """
def CalcChi1n( mol: Mol, force: bool = False) -> float:
    """
    CalcChi1n( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi1n(RDKit::ROMol [,bool=False])
    """
def CalcChi1v( mol: Mol, force: bool = False) -> float:
    """
    CalcChi1v( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi1v(RDKit::ROMol [,bool=False])
    """
def CalcChi2n( mol: Mol, force: bool = False) -> float:
    """
    CalcChi2n( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi2n(RDKit::ROMol [,bool=False])
    """
def CalcChi2v( mol: Mol, force: bool = False) -> float:
    """
    CalcChi2v( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi2v(RDKit::ROMol [,bool=False])
    """
def CalcChi3n( mol: Mol, force: bool = False) -> float:
    """
    CalcChi3n( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi3n(RDKit::ROMol [,bool=False])
    """
def CalcChi3v( mol: Mol, force: bool = False) -> float:
    """
    CalcChi3v( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi3v(RDKit::ROMol [,bool=False])
    """
def CalcChi4n( mol: Mol, force: bool = False) -> float:
    """
    CalcChi4n( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi4n(RDKit::ROMol [,bool=False])
    """
def CalcChi4v( mol: Mol, force: bool = False) -> float:
    """
    CalcChi4v( mol: Mol, force: bool = False) -> float

        C++ signature :
            double CalcChi4v(RDKit::ROMol [,bool=False])
    """
def CalcChiNn( mol: Mol, n: int, force: bool = False) -> float:
    """
    CalcChiNn( mol: Mol, n: int, force: bool = False) -> float

        C++ signature :
            double CalcChiNn(RDKit::ROMol,unsigned int [,bool=False])
    """
def CalcChiNv( mol: Mol, n: int, force: bool = False) -> float:
    """
    CalcChiNv( mol: Mol, n: int, force: bool = False) -> float

        C++ signature :
            double CalcChiNv(RDKit::ROMol,unsigned int [,bool=False])
    """
def CalcCoulombMat( mol: Mol, confId: int = -1) -> tuple:
    """
    CalcCoulombMat( mol: Mol, confId: int = -1) -> tuple
        Returns severals Coulomb randomized matrices

        C++ signature :
            boost::python::tuple CalcCoulombMat(RDKit::ROMol [,int=-1])
    """
def CalcCrippenDescriptors( mol: Mol, includeHs: bool = True, force: bool = False) -> tuple:
    """
    CalcCrippenDescriptors( mol: Mol, includeHs: bool = True, force: bool = False) -> tuple
        returns a 2-tuple with the Wildman-Crippen logp,mr values

        C++ signature :
            boost::python::tuple CalcCrippenDescriptors(RDKit::ROMol [,bool=True [,bool=False]])
    """
def CalcEEMcharges( mol: Mol, confId: int = -1) -> list:
    """
    CalcEEMcharges( mol: Mol, confId: int = -1) -> list
        Returns EEM atomic partial charges

        C++ signature :
            boost::python::list CalcEEMcharges(RDKit::ROMol {lvalue} [,int=-1])
    """
def CalcEccentricity( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcEccentricity( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcEccentricity(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcExactMolWt( mol: Mol, onlyHeavy: bool = False) -> float:
    """
    CalcExactMolWt( mol: Mol, onlyHeavy: bool = False) -> float
        returns the molecule's exact molecular weight

        C++ signature :
            double CalcExactMolWt(RDKit::ROMol [,bool=False])
    """
def CalcFractionCSP3( mol: Mol) -> float:
    """
    CalcFractionCSP3( mol: Mol) -> float
        returns the fraction of C atoms that are SP3 hybridized

        C++ signature :
            double CalcFractionCSP3(RDKit::ROMol)
    """
def CalcGETAWAY( mol: Mol, confId: int = -1, precision: float = 2, CustomAtomProperty: str = '') -> list:
    """
    CalcGETAWAY( mol: Mol, confId: int = -1, precision: float = 2, CustomAtomProperty: str = '') -> list
        Returns the GETAWAY descriptors vector

        C++ signature :
            boost::python::list CalcGETAWAY(RDKit::ROMol [,int=-1 [,double=2 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']]])
    """
def CalcHallKierAlpha( mol: Mol, atomContribs: AtomPairsParameters = None) -> float:
    """
    CalcHallKierAlpha( mol: Mol, atomContribs: AtomPairsParameters = None) -> float

        C++ signature :
            double CalcHallKierAlpha(RDKit::ROMol [,boost::python::api::object=None])
    """
def CalcInertialShapeFactor( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcInertialShapeFactor( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcInertialShapeFactor(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcKappa1( mol: Mol) -> float:
    """
    CalcKappa1( mol: Mol) -> float

        C++ signature :
            double CalcKappa1(RDKit::ROMol)
    """
def CalcKappa2( mol: Mol) -> float:
    """
    CalcKappa2( mol: Mol) -> float

        C++ signature :
            double CalcKappa2(RDKit::ROMol)
    """
def CalcKappa3( mol: Mol) -> float:
    """
    CalcKappa3( mol: Mol) -> float

        C++ signature :
            double CalcKappa3(RDKit::ROMol)
    """
def CalcLabuteASA( mol: Mol, includeHs: bool = True, force: bool = False) -> float:
    """
    CalcLabuteASA( mol: Mol, includeHs: bool = True, force: bool = False) -> float
        returns the Labute ASA value for a molecule

        C++ signature :
            double CalcLabuteASA(RDKit::ROMol [,bool=True [,bool=False]])
    """
def CalcMORSE( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
    CalcMORSE( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list
        Returns Molecule Representation of Structures based on Electron diffraction descriptors

        C++ signature :
            boost::python::list CalcMORSE(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']])
    """
def CalcMolFormula( mol: Mol, separateIsotopes: bool = False, abbreviateHIsotopes: bool = True) -> str:
    """
    CalcMolFormula( mol: Mol, separateIsotopes: bool = False, abbreviateHIsotopes: bool = True) -> str
        returns the molecule's formula

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > CalcMolFormula(RDKit::ROMol [,bool=False [,bool=True]])
    """
def CalcNPR1( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcNPR1( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcNPR1(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcNPR2( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcNPR2( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcNPR2(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcNumAliphaticCarbocycles( mol: Mol) -> int:
    """
    CalcNumAliphaticCarbocycles( mol: Mol) -> int
        returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticCarbocycles(RDKit::ROMol)
    """
def CalcNumAliphaticHeterocycles( mol: Mol) -> int:
    """
    CalcNumAliphaticHeterocycles( mol: Mol) -> int
        returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticHeterocycles(RDKit::ROMol)
    """
def CalcNumAliphaticRings( mol: Mol) -> int:
    """
    CalcNumAliphaticRings( mol: Mol) -> int
        returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticRings(RDKit::ROMol)
    """
def CalcNumAmideBonds( mol: Mol) -> int:
    """
    CalcNumAmideBonds( mol: Mol) -> int
        returns the number of amide bonds in a molecule

        C++ signature :
            unsigned int CalcNumAmideBonds(RDKit::ROMol)
    """
def CalcNumAromaticCarbocycles( mol: Mol) -> int:
    """
    CalcNumAromaticCarbocycles( mol: Mol) -> int
        returns the number of aromatic carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumAromaticCarbocycles(RDKit::ROMol)
    """
def CalcNumAromaticHeterocycles( mol: Mol) -> int:
    """
    CalcNumAromaticHeterocycles( mol: Mol) -> int
        returns the number of aromatic heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumAromaticHeterocycles(RDKit::ROMol)
    """
def CalcNumAromaticRings( mol: Mol) -> int:
    """
    CalcNumAromaticRings( mol: Mol) -> int
        returns the number of aromatic rings for a molecule

        C++ signature :
            unsigned int CalcNumAromaticRings(RDKit::ROMol)
    """
def CalcNumAtomStereoCenters( mol: Mol) -> int:
    """
    CalcNumAtomStereoCenters( mol: Mol) -> int
        Returns the total number of atomic stereocenters (specified and unspecified)

        C++ signature :
            unsigned int CalcNumAtomStereoCenters(RDKit::ROMol)
    """
def CalcNumAtoms( mol: Mol) -> int:
    """
    CalcNumAtoms( mol: Mol) -> int
        returns the total number of atoms for a molecule

        C++ signature :
            unsigned int CalcNumAtoms(RDKit::ROMol)
    """
def CalcNumBridgeheadAtoms( mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
    CalcNumBridgeheadAtoms( mol: Mol, atoms: AtomPairsParameters = None) -> int
        Returns the number of bridgehead atoms (atoms shared between rings that share at least two bonds)

        C++ signature :
            unsigned int CalcNumBridgeheadAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
def CalcNumHBA( mol: Mol) -> int:
    """
    CalcNumHBA( mol: Mol) -> int
        returns the number of H-bond acceptors for a molecule

        C++ signature :
            unsigned int CalcNumHBA(RDKit::ROMol)
    """
def CalcNumHBD( mol: Mol) -> int:
    """
    CalcNumHBD( mol: Mol) -> int
        returns the number of H-bond donors for a molecule

        C++ signature :
            unsigned int CalcNumHBD(RDKit::ROMol)
    """
def CalcNumHeavyAtoms( mol: Mol) -> int:
    """
    CalcNumHeavyAtoms( mol: Mol) -> int
        returns the number of heavy atoms for a molecule

        C++ signature :
            unsigned int CalcNumHeavyAtoms(RDKit::ROMol)
    """
def CalcNumHeteroatoms( mol: Mol) -> int:
    """
    CalcNumHeteroatoms( mol: Mol) -> int
        returns the number of heteroatoms for a molecule

        C++ signature :
            unsigned int CalcNumHeteroatoms(RDKit::ROMol)
    """
def CalcNumHeterocycles( mol: Mol) -> int:
    """
    CalcNumHeterocycles( mol: Mol) -> int
        returns the number of heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumHeterocycles(RDKit::ROMol)
    """
def CalcNumLipinskiHBA( mol: Mol) -> int:
    """
    CalcNumLipinskiHBA( mol: Mol) -> int
        returns the number of Lipinski H-bond acceptors for a molecule

        C++ signature :
            unsigned int CalcNumLipinskiHBA(RDKit::ROMol)
    """
def CalcNumLipinskiHBD( mol: Mol) -> int:
    """
    CalcNumLipinskiHBD( mol: Mol) -> int
        returns the number of Lipinski H-bond donors for a molecule

        C++ signature :
            unsigned int CalcNumLipinskiHBD(RDKit::ROMol)
    """
def CalcNumRings( mol: Mol) -> int:
    """
    CalcNumRings( mol: Mol) -> int
        returns the number of rings for a molecule

        C++ signature :
            unsigned int CalcNumRings(RDKit::ROMol)
    """
@typing.overload
def CalcNumRotatableBonds( mol: Mol, strict: bool) -> int:
    """
    CalcNumRotatableBonds( mol: Mol, strict: bool) -> int
        returns the number of rotatable bonds for a molecule.
           strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
           strict = NumRotatableBondsOptions.Strict - (default) does not count things like
                    amide or ester bonds
           strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
              systems.
              - Single bonds between aliphatic ring Cs are always rotatable. This
                means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now 
                considered rotatable; it was not before
              - Heteroatoms in the linked rings no longer affect whether or not
                the linking bond is rotatable
              - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
                 considered non-rotatable

        C++ signature :
            unsigned int CalcNumRotatableBonds(RDKit::ROMol,bool)

        C++ signature :
            unsigned int CalcNumRotatableBonds(RDKit::ROMol [,RDKit::Descriptors::NumRotatableBondsOptions=rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default])
    """
@typing.overload
def CalcNumRotatableBonds( mol: Mol, strict: NumRotatableBondsOptions = rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default) -> int:
    pass
def CalcNumSaturatedCarbocycles( mol: Mol) -> int:
    """
    CalcNumSaturatedCarbocycles( mol: Mol) -> int
        returns the number of saturated carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedCarbocycles(RDKit::ROMol)
    """
def CalcNumSaturatedHeterocycles( mol: Mol) -> int:
    """
    CalcNumSaturatedHeterocycles( mol: Mol) -> int
        returns the number of saturated heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedHeterocycles(RDKit::ROMol)
    """
def CalcNumSaturatedRings( mol: Mol) -> int:
    """
    CalcNumSaturatedRings( mol: Mol) -> int
        returns the number of saturated rings for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedRings(RDKit::ROMol)
    """
def CalcNumSpiroAtoms( mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
    CalcNumSpiroAtoms( mol: Mol, atoms: AtomPairsParameters = None) -> int
        Returns the number of spiro atoms (atoms shared between rings that share exactly one atom)

        C++ signature :
            unsigned int CalcNumSpiroAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
def CalcNumUnspecifiedAtomStereoCenters( mol: Mol) -> int:
    """
    CalcNumUnspecifiedAtomStereoCenters( mol: Mol) -> int
        Returns the number of unspecified atomic stereocenters

        C++ signature :
            unsigned int CalcNumUnspecifiedAtomStereoCenters(RDKit::ROMol)
    """
def CalcOxidationNumbers( mol: Mol) -> None:
    """
    CalcOxidationNumbers( mol: Mol) -> None
        Adds the oxidation number/state to the atoms of a molecule as property OxidationNumber on each atom.  Use Pauling electronegativities.  This is experimental code, still under development.

        C++ signature :
            void CalcOxidationNumbers(RDKit::ROMol)
    """
def CalcPBF( mol: Mol, confId: int = -1) -> float:
    """
    CalcPBF( mol: Mol, confId: int = -1) -> float
        Returns the PBF (plane of best fit) descriptor (https://doi.org/10.1021/ci300293f)

        C++ signature :
            double CalcPBF(RDKit::ROMol [,int=-1])
    """
def CalcPMI1( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcPMI1( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcPMI1(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcPMI2( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcPMI2( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcPMI2(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcPMI3( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcPMI3( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcPMI3(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcPhi( mol: Mol) -> float:
    """
    CalcPhi( mol: Mol) -> float

        C++ signature :
            double CalcPhi(RDKit::ROMol)
    """
def CalcRDF( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
    CalcRDF( mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list
        Returns radial distribution fonction descriptors (RDF)

        C++ signature :
            boost::python::list CalcRDF(RDKit::ROMol [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']])
    """
def CalcRMS( prbMol: Mol, refMol: Mol, prbId: int = -1, refId: int = -1, map: AtomPairsParameters = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = []) -> float:
    """
    CalcRMS( prbMol: Mol, refMol: Mol, prbId: int = -1, refId: int = -1, map: AtomPairsParameters = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = []) -> float
        Returns the RMS between two molecules, taking symmetry into account.
               In contrast to getBestRMS, the RMS is computed 'in place', i.e.
               probe molecules are not aligned to the reference ahead of the
               RMS calculation. This is useful, for example, to compute
               the RMSD between docking poses and the co-crystallized ligand.
              
               Note:
               This function will attempt to match all permutations of matching atom
               orders in both molecules, for some molecules it will lead to
               'combinatorial explosion' especially if hydrogens are present.
              
               ARGUMENTS
                - prbMol:      the molecule to be aligned to the reference
                - refMol:      the reference molecule
                - prbCId:      (optional) probe conformation to use
                - refCId:      (optional) reference conformation to use
                - map:         (optional) a list of lists of (probeAtomId, refAtomId)
                               tuples with the atom-atom mappings of the two
                               molecules. If not provided, these will be generated
                               using a substructure search.
                - maxMatches:  (optional) if map isn't specified, this will be
                               the max number of matches found in a SubstructMatch()
                - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
                               terminal functional groups (like nitro or carboxylate)
                               will be considered symmetrically
                - weights:     (optional) weights for mapping
               
              RETURNS
              The best RMSD found
            
        

        C++ signature :
            double CalcRMS(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1 [,boost::python::api::object=None [,int=1000000 [,bool=True [,boost::python::api::object=[]]]]]]])
    """
def CalcRadiusOfGyration( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    CalcRadiusOfGyration( mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float

        C++ signature :
            double CalcRadiusOfGyration(RDKit::ROMol [,int=-1 [,bool=True [,bool=True]]])
    """
def CalcSpherocityIndex( mol: Mol, confId: int = -1, force: bool = True) -> float:
    """
    CalcSpherocityIndex( mol: Mol, confId: int = -1, force: bool = True) -> float

        C++ signature :
            double CalcSpherocityIndex(RDKit::ROMol [,int=-1 [,bool=True]])
    """
def CalcTPSA( mol: Mol, force: bool = False, includeSandP: bool = False) -> float:
    """
    CalcTPSA( mol: Mol, force: bool = False, includeSandP: bool = False) -> float
        returns the TPSA value for a molecule

        C++ signature :
            double CalcTPSA(RDKit::ROMol [,bool=False [,bool=False]])
    """
def CalcWHIM( mol: Mol, confId: int = -1, thresh: float = 0.001, CustomAtomProperty: str = '') -> list:
    """
    CalcWHIM( mol: Mol, confId: int = -1, thresh: float = 0.001, CustomAtomProperty: str = '') -> list
        Returns the WHIM descriptors vector

        C++ signature :
            boost::python::list CalcWHIM(RDKit::ROMol [,int=-1 [,double=0.001 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']]])
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
def CanonicalRankAtomsInFragment( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vecti:
    """
    CanonicalRankAtomsInFragment( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, breakTies: bool = True, includeChirality: bool = True, includeIsotopes: bool = True) -> _vecti
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
def CanonicalizeConformer( conf: Conformer, center: Point3D = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> None:
    """
    CanonicalizeConformer( conf: Conformer, center: Point3D = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> None
        Canonicalize the orientation of a conformer so that its principal axes
                       around the specified center point coincide with the x, y, z axes
          
          ARGUMENTS:
            - conf : conformer of interest 
            - center : optionally center point about which the principal axes are computed 
                                  if not specified the centroid of the conformer will be used
            - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms
        

        C++ signature :
            void CanonicalizeConformer(RDKit::Conformer {lvalue} [,RDGeom::Point3D const*=None [,bool=False [,bool=True]]])
    """
def CanonicalizeEnhancedStereo( mol: Mol) -> None:
    """
    CanonicalizeEnhancedStereo( mol: Mol) -> None

        C++ signature :
            void CanonicalizeEnhancedStereo(RDKit::ROMol {lvalue})
    """
def CanonicalizeMol( mol: Mol, normalizeCovar: bool = False, ignoreHs: bool = True) -> None:
    """
    CanonicalizeMol( mol: Mol, normalizeCovar: bool = False, ignoreHs: bool = True) -> None
        Loop over the conformers in a molecule and canonicalize their orientation

        C++ signature :
            void CanonicalizeMol(RDKit::ROMol {lvalue} [,bool=False [,bool=True]])
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
def Compute2DCoords( mol: Mol, canonOrient: bool = True, clearConfs: bool = True, coordMap: dict = {}, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0, forceRDKit: bool = False, useRingTemplates: bool = False) -> int:
    """
    Compute2DCoords( mol: Mol, canonOrient: bool = True, clearConfs: bool = True, coordMap: dict = {}, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0, forceRDKit: bool = False, useRingTemplates: bool = False) -> int
        Compute 2D coordinates for a molecule. 
          The resulting coordinates are stored on each atom of the molecule 
        
          ARGUMENTS: 
        
             mol - the molecule of interest
             canonOrient - orient the molecule in a canonical way
             clearConfs - if true, all existing conformations on the molecule
                     will be cleared
             coordMap - a dictionary mapping atom Ids -> Point2D objects 
                        with starting coordinates for atoms that should
                        have their positions locked.
             nFlipsPerSample - number of rotatable bonds that are
                        flipped at random at a time.
             nSample - Number of random samplings of rotatable bonds.
             sampleSeed - seed for the random sampling process.
             permuteDeg4Nodes - allow permutation of bonds at a degree 4
                         node during the sampling process 
             bondLength - change the default bond length for depiction 
             forceRDKit - use RDKit to generate coordinates even if 
                          preferCoordGen is set to true
             useRingTemplates - use templates to generate coordinates of complex
                          ring systems
        
          RETURNS: 
        
             ID of the conformation added to the molecule
        

        C++ signature :
            unsigned int Compute2DCoords(RDKit::ROMol {lvalue} [,bool=True [,bool=True [,boost::python::dict {lvalue}={} [,unsigned int=0 [,unsigned int=0 [,int=0 [,bool=False [,double=-1.0 [,bool=False [,bool=False]]]]]]]]]])
    """
def Compute2DCoordsForReaction( reaction: ChemicalReaction, spacing: float = 2.0, updateProps: bool = True, canonOrient: bool = True, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0) -> None:
    """
    Compute2DCoordsForReaction( reaction: ChemicalReaction, spacing: float = 2.0, updateProps: bool = True, canonOrient: bool = True, nFlipsPerSample: int = 0, nSample: int = 0, sampleSeed: int = 0, permuteDeg4Nodes: bool = False, bondLength: float = -1.0) -> None
        Compute 2D coordinates for a reaction. 
          ARGUMENTS: 
             - reaction - the reaction of interest
             - spacing - the amount of space left between components of the reaction
             - canonOrient - orient the reactants and products in a canonical way
             - updateProps - if set, properties such as conjugation and
                hybridization will be calculated for the reactant and product
                templates before generating coordinates. This should result in
                better depictions, but can lead to errors in some cases.
             - nFlipsPerSample - number of rotatable bonds that are
                        flipped at random at a time.
             - nSample - Number of random samplings of rotatable bonds.
             - sampleSeed - seed for the random sampling process.
             - permuteDeg4Nodes - allow permutation of bonds at a degree 4
                         node during the sampling process 
             - bondLength - change the default bond length for depiction
        

        C++ signature :
            void Compute2DCoordsForReaction(RDKit::ChemicalReaction {lvalue} [,double=2.0 [,bool=True [,bool=True [,unsigned int=0 [,unsigned int=0 [,int=0 [,bool=False [,double=-1.0]]]]]]]])
    """
def Compute2DCoordsMimicDistmat( mol: Mol, distMat: AtomPairsParameters, canonOrient: bool = False, clearConfs: bool = True, weightDistMat: float = 0.5, nFlipsPerSample: int = 3, nSample: int = 100, sampleSeed: int = 100, permuteDeg4Nodes: bool = True, bondLength: float = -1.0, forceRDKit: bool = False) -> int:
    """
    Compute2DCoordsMimicDistmat( mol: Mol, distMat: AtomPairsParameters, canonOrient: bool = False, clearConfs: bool = True, weightDistMat: float = 0.5, nFlipsPerSample: int = 3, nSample: int = 100, sampleSeed: int = 100, permuteDeg4Nodes: bool = True, bondLength: float = -1.0, forceRDKit: bool = False) -> int
        Compute 2D coordinates for a molecule such 
          that the inter-atom distances mimic those in a user-provided
          distance matrix. 
          The resulting coordinates are stored on each atom of the molecule 
        
          ARGUMENTS: 
        
             mol - the molecule of interest
             distMat - distance matrix that we want the 2D structure to mimic
             canonOrient - orient the molecule in a canonical way
             clearConfs - if true, all existing conformations on the molecule
                     will be cleared
             weightDistMat - weight assigned in the cost function to mimicking
                             the distance matrix.
                             This must be between (0.0,1.0). (1.0-weightDistMat)
                             is then the weight assigned to improving 
                             the density of the 2D structure i.e. try to
                             make it spread out
             nFlipsPerSample - number of rotatable bonds that are
                        flipped at random at a time.
             nSample - Number of random samplings of rotatable bonds.
             sampleSeed - seed for the random sampling process.
             permuteDeg4Nodes - allow permutation of bonds at a degree 4
                         node during the sampling process 
             bondLength - change the default bond length for depiction 
             forceRDKit - use RDKit to generate coordinates even if 
                          preferCoordGen is set to true
        
          RETURNS: 
        
             ID of the conformation added to the molecule
        

        C++ signature :
            unsigned int Compute2DCoordsMimicDistmat(RDKit::ROMol {lvalue},boost::python::api::object [,bool=False [,bool=True [,double=0.5 [,unsigned int=3 [,unsigned int=100 [,int=100 [,bool=True [,double=-1.0 [,bool=False]]]]]]]]])
    """
def ComputeCanonicalTransform( conf: Conformer, center: Point3D = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> object:
    """
    ComputeCanonicalTransform( conf: Conformer, center: Point3D = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> object
        Compute the transformation required aligna conformer so that
                       the principal axes align up with the x,y, z axes
                       The conformer itself is left unchanged
          ARGUMENTS:
            - conf : the conformer of interest
            - center : optional center point to compute the principal axes around (defaults to the centroid)
            - normalizeCovar : optionally normalize the covariance matrix by the number of atoms
        

        C++ signature :
            _object* ComputeCanonicalTransform(RDKit::Conformer [,RDGeom::Point3D const*=None [,bool=False [,bool=True]]])
    """
def ComputeCentroid( conf: Conformer, ignoreHs: bool = True, weights: _vectd = None) -> Point3D:
    """
    ComputeCentroid( conf: Conformer, ignoreHs: bool = True, weights: _vectd = None) -> Point3D
        Compute the centroid of the conformation - hydrogens are ignored and no attention
                                   is paid to the difference in sizes of the heavy atoms; however,
                                   an optional vector of weights can be passed.
        

        C++ signature :
            RDGeom::Point3D ComputeCentroid(RDKit::Conformer [,bool=True [,std::vector<double, std::allocator<double> > const*=None]])
    """
def ComputeConfBox( conf: Conformer, trans: AtomPairsParameters = None, padding: float = 2.0) -> tuple:
    """
    ComputeConfBox( conf: Conformer, trans: AtomPairsParameters = None, padding: float = 2.0) -> tuple
        Compute the lower and upper corners of a cuboid that will fit the conformer

        C++ signature :
            boost::python::tuple ComputeConfBox(RDKit::Conformer [,boost::python::api::object=None [,double=2.0]])
    """
def ComputeConfDimsAndOffset( conf: Conformer, trans: AtomPairsParameters = None, padding: float = 2.0) -> tuple:
    """
    ComputeConfDimsAndOffset( conf: Conformer, trans: AtomPairsParameters = None, padding: float = 2.0) -> tuple
        Compute the size of the box that can fit the conformations, and offset 
           of the box from the origin
        

        C++ signature :
            boost::python::tuple ComputeConfDimsAndOffset(RDKit::Conformer [,boost::python::api::object=None [,double=2.0]])
    """
def ComputeGasteigerCharges( mol: Mol, nIter: int = 12, throwOnParamFailure: bool = False) -> None:
    """
    ComputeGasteigerCharges( mol: Mol, nIter: int = 12, throwOnParamFailure: bool = False) -> None
        Compute Gasteiger partial charges for molecule
        
         The charges are computed using an iterative procedure presented in 
         
         Ref : J.Gasteiger, M. Marseli, Iterative Equalization of Oribital Electronegatiity 
         A Rapid Access to Atomic Charges, Tetrahedron Vol 36 p3219 1980
         
         The computed charges are stored on each atom are stored a computed property ( under the name 
         _GasteigerCharge). In addition, each atom also stored the total charge for the implicit hydrogens 
         on the atom (under the property name _GasteigerHCharge)
         
         ARGUMENTS:
        
            - mol : the molecule of interrest
            - nIter : number of iteration (defaults to 12)
            - throwOnParamFailure : toggles whether or not an exception should be raised if parameters
              for an atom cannot be found.  If this is false (the default), all parameters for unknown
              atoms will be set to zero.  This has the effect of removing that atom from the iteration.
        
        

        C++ signature :
            void ComputeGasteigerCharges(RDKit::ROMol [,int=12 [,bool=False]])
    """
def ComputePrincipalAxesAndMoments( conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None) -> object:
    """
    ComputePrincipalAxesAndMoments( conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None) -> object
        Compute principal axes and moments of inertia for a conformer
               These values are calculated from the inertia tensor:
               Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
               Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})
               where the coordinates are relative to the center of mass.
        
          ARGUMENTS:
            - conf : the conformer of interest
            - ignoreHs : if True, ignore hydrogen atoms
            - weights : if present, used to weight the atomic coordinates
        
          Returns a (principal axes, principal moments) tuple
        

        C++ signature :
            _object* ComputePrincipalAxesAndMoments(RDKit::Conformer [,bool=True [,boost::python::api::object=None]])
    """
def ComputePrincipalAxesAndMomentsFromGyrationMatrix( conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None) -> object:
    """
    ComputePrincipalAxesAndMomentsFromGyrationMatrix( conf: Conformer, ignoreHs: bool = True, weights: AtomPairsParameters = None) -> object
        Compute principal axes and moments from the gyration matrix of a conformer
               These values are calculated from the gyration matrix/tensor:
               Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
               Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})
               where the coordinates are relative to the center of mass.
        
          ARGUMENTS:
            - conf : the conformer of interest
            - ignoreHs : if True, ignore hydrogen atoms
            - weights : if present, used to weight the atomic coordinates
        
          Returns a (principal axes, principal moments) tuple
        

        C++ signature :
            _object* ComputePrincipalAxesAndMomentsFromGyrationMatrix(RDKit::Conformer [,bool=True [,boost::python::api::object=None]])
    """
def ComputeUnionBox( arg1: tuple, arg2: tuple) -> tuple:
    """
    ComputeUnionBox( arg1: tuple, arg2: tuple) -> tuple
        Compute the union of two boxes, so that all the points in both boxes are 
            contained in the new box

        C++ signature :
            boost::python::tuple ComputeUnionBox(boost::python::tuple,boost::python::tuple)
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
def CreateDifferenceFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> UIntSparseIntVect:
    """
    CreateDifferenceFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> UIntSparseIntVect
        construct a difference fingerprint for a ChemicalReaction by subtracting the reactant fingerprint from the product fingerprint

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* CreateDifferenceFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7f2c64d21840>])
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
def CreateStereoGroup( stereoGroupType: StereoGroupType, mol: Mol, atomIds: AtomPairsParameters, readId: int = 0) -> StereoGroup:
    """
    CreateStereoGroup( stereoGroupType: StereoGroupType, mol: Mol, atomIds: AtomPairsParameters, readId: int = 0) -> StereoGroup
        creates a StereoGroup associated with a molecule from a list of atom Ids

        C++ signature :
            RDKit::StereoGroup* CreateStereoGroup(RDKit::StereoGroupType,RDKit::ROMol {lvalue},boost::python::api::object [,unsigned int=0])
    """
def CreateStructuralFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> ExplicitBitVect:
    """
    CreateStructuralFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> ExplicitBitVect
        construct a structural fingerprint for a ChemicalReaction by concatenating the reactant fingerprint and the product fingerprint

        C++ signature :
            ExplicitBitVect* CreateStructuralFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7f2c64d21240>])
    """
def CustomProp_VSA_( mol: Mol, customPropName: str, bins: AtomPairsParameters, force: bool = False) -> list:
    """
    CustomProp_VSA_( mol: Mol, customPropName: str, bins: AtomPairsParameters, force: bool = False) -> list

        C++ signature :
            boost::python::list CustomProp_VSA_(RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,boost::python::api::object [,bool=False])
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
def DetectChemistryProblems( mol: Mol, sanitizeOps: int = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL) -> tuple:
    """
    DetectChemistryProblems( mol: Mol, sanitizeOps: int = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL) -> tuple
        checks for chemistry problems

        C++ signature :
            boost::python::tuple DetectChemistryProblems(RDKit::ROMol [,unsigned int=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL])
    """
def ETDG() -> EmbedParameters:
    """
    ETDG() -> EmbedParameters
        Returns an EmbedParameters object for the ETDG method.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETDG()
    """
def ETKDG() -> EmbedParameters:
    """
    ETKDG() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 1.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDG()
    """
def ETKDGv2() -> EmbedParameters:
    """
    ETKDGv2() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 2.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDGv2()
    """
def ETKDGv3() -> EmbedParameters:
    """
    ETKDGv3() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 3 (macrocycles).

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDGv3()
    """
@typing.overload
def EmbedMolecule( mol: Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> int:
    """
    EmbedMolecule( mol: Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> int
        Use distance geometry to obtain initial 
         coordinates for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - maxAttempts : the maximum number of attempts to try embedding 
            - randomSeed : provide a seed for the random number generator 
                           so that the same coordinates can be obtained 
                           for a molecule on multiple runs. If -1, the 
                           RNG will not be seeded. 
            - clearConfs : clear all existing conformations on the molecule
            - useRandomCoords : Start the embedding from random coordinates instead of
                                using eigenvalues of the distance matrix.
            - boxSizeMult    Determines the size of the box that is used for
                             random coordinates. If this is a positive number, the 
                             side length will equal the largest element of the distance
                             matrix times boxSizeMult. If this is a negative number,
                             the side length will equal -boxSizeMult (i.e. independent
                             of the elements of the distance matrix).
            - randNegEig : If the embedding yields a negative eigenvalue, 
                           pick coordinates that correspond 
                           to this component at random 
            - numZeroFail : fail embedding if we have at least this many zero eigenvalues 
            - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                         require some atoms to have fixed coordinates in the resulting 
                         conformation.
            - forceTol : tolerance to be used during the force-field minimization with 
                         the distance geometry force field.
            - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                         of the bounds matrix fails.
            - enforceChirality : enforce the correct chirality if chiral centers are present.
            - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
            - useBasicKnowledge : impose basic knowledge such as flat rings
            - printExpTorsionAngles : print the output from the experimental torsion angles
        
         RETURNS:
        
            ID of the new conformation added to the molecule 
        
        

        C++ signature :
            int EmbedMolecule(RDKit::ROMol {lvalue} [,unsigned int=0 [,int=-1 [,bool=True [,bool=False [,double=2.0 [,bool=True [,unsigned int=1 [,boost::python::dict {lvalue}={} [,double=0.001 [,bool=False [,bool=True [,bool=True [,bool=True [,bool=False [,bool=False [,bool=False [,unsigned int=1]]]]]]]]]]]]]]]]])

        C++ signature :
            int EmbedMolecule(RDKit::ROMol {lvalue},RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
@typing.overload
def EmbedMolecule( mol: Mol, params: EmbedParameters) -> int:
    pass
@typing.overload
def EmbedMultipleConfs( mol: Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> _vecti:
    """
    EmbedMultipleConfs( mol: Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> _vecti
        Use distance geometry to obtain multiple sets of 
         coordinates for a molecule
         
         ARGUMENTS:
        
          - mol : the molecule of interest
          - numConfs : the number of conformers to generate 
          - maxAttempts : the maximum number of attempts to try embedding 
          - randomSeed : provide a seed for the random number generator 
                         so that the same coordinates can be obtained 
                         for a molecule on multiple runs. If -1, the 
                         RNG will not be seeded. 
          - clearConfs : clear all existing conformations on the molecule
          - useRandomCoords : Start the embedding from random coordinates instead of
                              using eigenvalues of the distance matrix.
          - boxSizeMult    Determines the size of the box that is used for
                           random coordinates. If this is a positive number, the 
                           side length will equal the largest element of the distance
                           matrix times boxSizeMult. If this is a negative number,
                           the side length will equal -boxSizeMult (i.e. independent
                           of the elements of the distance matrix).
          - randNegEig : If the embedding yields a negative eigenvalue, 
                         pick coordinates that correspond 
                         to this component at random 
          - numZeroFail : fail embedding if we have at least this many zero eigenvalues 
          - pruneRmsThresh : Retain only the conformations out of 'numConfs' 
                            after embedding that are at least 
                            this far apart from each other. 
                            RMSD is computed on the heavy atoms. 
                            Pruning is greedy; i.e. the first embedded conformation
                            is retained and from then on only those that are at
                            least pruneRmsThresh away from all retained conformations
                            are kept. The pruning is done after embedding and 
                            bounds violation minimization. No pruning by default.
          - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                       require some atoms to have fixed coordinates in the resulting 
                       conformation.
          - forceTol : tolerance to be used during the force-field minimization with 
                       the distance geometry force field.
          - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                       of the bounds matrix fails.
          - enforceChirality : enforce the correct chirality if chiral centers are present.
          - numThreads : number of threads to use while embedding. This only has an effect if the RDKit
                       was built with multi-thread support.
                      If set to zero, the max supported by the system will be used.
          - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
          - useBasicKnowledge : impose basic knowledge such as flat rings
          - printExpTorsionAngles : print the output from the experimental torsion angles
         RETURNS:
        
            List of new conformation IDs 
        
        

        C++ signature :
            std::vector<int, std::allocator<int> > EmbedMultipleConfs(RDKit::ROMol {lvalue} [,unsigned int=10 [,unsigned int=0 [,int=-1 [,bool=True [,bool=False [,double=2.0 [,bool=True [,unsigned int=1 [,double=-1.0 [,boost::python::dict {lvalue}={} [,double=0.001 [,bool=False [,bool=True [,int=1 [,bool=True [,bool=True [,bool=False [,bool=False [,bool=False [,unsigned int=1]]]]]]]]]]]]]]]]]]]])

        C++ signature :
            std::vector<int, std::allocator<int> > EmbedMultipleConfs(RDKit::ROMol {lvalue},unsigned int,RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
@typing.overload
def EmbedMultipleConfs( mol: Mol, numConfs: int, params: EmbedParameters) -> _vecti:
    pass
def EncodeShape( mol: Mol, grid: UniformGrid3D_, confId: int = -1, trans: AtomPairsParameters = None, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> None:
    """
    EncodeShape( mol: Mol, grid: UniformGrid3D_, confId: int = -1, trans: AtomPairsParameters = None, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> None
        Encode the shape of a molecule (one of its conformer) onto a grid
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - grid : grid onto which the encoding is written 
            - confId : id of the conformation of interest on mol (defaults to the first one) 
            - trans : any transformation that needs to be used to encode onto the grid (note the molecule remains unchanged) 
            - vdwScale : Scaling factor for the radius of the atoms to determine the base radius 
                         used in the encoding - grid points inside this sphere carry the maximum occupancy 
            - setpSize : thickness of the layers outside the base radius, the occupancy value is decreased 
                         from layer to layer from the maximum value 
            - maxLayers : the maximum number of layers - defaults to the number of bits 
                          used per grid point - e.g. two bits per grid point will allow 3 layers
            - ignoreHs : when set, the contribution of Hs to the shape will be ignored
        

        C++ signature :
            void EncodeShape(RDKit::ROMol,RDGeom::UniformGrid3D {lvalue} [,int=-1 [,boost::python::api::object=None [,double=0.8 [,double=0.25 [,int=-1 [,bool=True]]]]]])
    """
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
def EnumerateLibraryCanSerialize() -> bool:
    """
    EnumerateLibraryCanSerialize() -> bool
        Returns True if the EnumerateLibrary is serializable (requires boost serialization

        C++ signature :
            bool EnumerateLibraryCanSerialize()
    """
def ExplicitDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeEqualsQueryAtom(int [,bool=False])
    """
def ExplicitDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeGreaterQueryAtom(int [,bool=False])
    """
def ExplicitDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeLessQueryAtom(int [,bool=False])
    """
def ExplicitValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceEqualsQueryAtom(int [,bool=False])
    """
def ExplicitValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceGreaterQueryAtom(int [,bool=False])
    """
def ExplicitValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceLessQueryAtom(int [,bool=False])
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
def FindAtomEnvironmentOfRadiusN( mol: Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: AtomPairsParameters = None) -> _vecti:
    """
    FindAtomEnvironmentOfRadiusN( mol: Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: AtomPairsParameters = None) -> _vecti
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
def FormalChargeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* FormalChargeEqualsQueryAtom(int [,bool=False])
    """
def FormalChargeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* FormalChargeGreaterQueryAtom(int [,bool=False])
    """
def FormalChargeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* FormalChargeLessQueryAtom(int [,bool=False])
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
def FragmentOnBonds( mol: Mol, bondIndices: AtomPairsParameters, addDummies: bool = True, dummyLabels: AtomPairsParameters = None, bondTypes: AtomPairsParameters = None, cutsPerAtom: list = []) -> Mol:
    """
    FragmentOnBonds( mol: Mol, bondIndices: AtomPairsParameters, addDummies: bool = True, dummyLabels: AtomPairsParameters = None, bondTypes: AtomPairsParameters = None, cutsPerAtom: list = []) -> Mol
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
def FragmentOnSomeBonds( mol: Mol, bondIndices: AtomPairsParameters, numToBreak: int = 1, addDummies: bool = True, dummyLabels: AtomPairsParameters = None, bondTypes: AtomPairsParameters = None, returnCutsPerAtom: bool = False) -> tuple:
    """
    FragmentOnSomeBonds( mol: Mol, bondIndices: AtomPairsParameters, numToBreak: int = 1, addDummies: bool = True, dummyLabels: AtomPairsParameters = None, bondTypes: AtomPairsParameters = None, returnCutsPerAtom: bool = False) -> tuple
        fragment on some bonds

        C++ signature :
            boost::python::tuple FragmentOnSomeBonds(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None [,bool=False]]]]])
    """
@typing.overload
def GenerateDepictionMatching2DStructure( mol: Mol, reference: Mol, confId: int = -1, refPatt: AtomPairsParameters = None, acceptFailure: bool = False, forceRDKit: bool = False, allowRGroups: bool = False) -> tuple:
    """
    GenerateDepictionMatching2DStructure( mol: Mol, reference: Mol, confId: int = -1, refPatt: AtomPairsParameters = None, acceptFailure: bool = False, forceRDKit: bool = False, allowRGroups: bool = False) -> tuple
        Generate a depiction for a molecule where a piece of the 
          molecule is constrained to have the same coordinates as a reference. 
        
          This is useful for, for example, generating depictions of SAR data 
          sets so that the cores of the molecules are all oriented the same way. 
          ARGUMENTS: 
        
          mol -    the molecule to be aligned, this will come back 
                   with a single conformer. 
          reference -    a molecule with the reference atoms to align to; 
                         this should have a depiction. 
          confId -       (optional) the id of the reference conformation to use 
          refPatt -      (optional) a query molecule to be used to generate 
                         the atom mapping between the molecule and the reference 
          acceptFailure - (optional) if True, standard depictions will be generated 
                          for molecules that don't have a substructure match to the 
                          reference; if False, throws a DepictException.
          forceRDKit -    (optional) use RDKit to generate coordinates even if 
                          preferCoordGen is set to true
          allowRGroups -  (optional) if True, terminal dummy atoms in the 
                          reference are ignored if they match an implicit 
                          hydrogen in the molecule, and a constrained 
                          depiction is still attempted
        
          RETURNS: a tuple of (refIdx, molIdx) tuples corresponding to the atom 
                   indices in mol constrained to have the same coordinates as atom 
                   indices in reference.
        

        C++ signature :
            boost::python::tuple GenerateDepictionMatching2DStructure(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,boost::python::api::object=None [,bool=False [,bool=False [,bool=False]]]]])

        C++ signature :
            void GenerateDepictionMatching2DStructure(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue},boost::python::api::object [,int=-1 [,bool=False]])
    """
@typing.overload
def GenerateDepictionMatching2DStructure( mol: Mol, reference: Mol, atomMap: AtomPairsParameters, confId: int = -1, forceRDKit: bool = False) -> None:
    pass
def GenerateDepictionMatching3DStructure( mol: Mol, reference: Mol, confId: int = -1, refPatt: AtomPairsParameters = None, acceptFailure: bool = False, forceRDKit: bool = False) -> None:
    """
    GenerateDepictionMatching3DStructure( mol: Mol, reference: Mol, confId: int = -1, refPatt: AtomPairsParameters = None, acceptFailure: bool = False, forceRDKit: bool = False) -> None
        Generate a depiction for a molecule where a piece of the molecule 
          is constrained to have coordinates similar to those of a 3D reference 
          structure.
          ARGUMENTS: 
        
          mol -    the molecule to be aligned, this will come back 
                   with a single conformer containing the 2D coordinates. 
          reference -    a molecule with the reference atoms to align to. 
                         By default this should be the same as mol, but with 
                         3D coordinates 
          confId -       (optional) the id of the reference conformation to use 
          referencePattern -  (optional) a query molecule to map a subset of 
                              the reference onto the mol, so that only some of the 
                              atoms are aligned. 
          acceptFailure - (optional) if True, standard depictions will be generated 
                          for molecules that don't match the reference or the
                          referencePattern; if False, throws a DepictException.
          forceRDKit -    (optional) use RDKit to generate coordinates even if 
                          preferCoordGen is set to true

        C++ signature :
            void GenerateDepictionMatching3DStructure(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,boost::python::api::object=None [,bool=False [,bool=False]]]])
    """
def GenerateErGFingerprintForReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object:
    """
    GenerateErGFingerprintForReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object
        Returns the ErG fingerprint vector for a reduced graph

        C++ signature :
            _object* GenerateErGFingerprintForReducedGraph(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
def GenerateMolExtendedReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0) -> Mol:
    """
    GenerateMolExtendedReducedGraph( mol: Mol, atomTypes: AtomPairsParameters = 0) -> Mol
        Returns the reduced graph for a molecule

        C++ signature :
            RDKit::ROMol* GenerateMolExtendedReducedGraph(RDKit::ROMol [,boost::python::api::object=0])
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
def GetAlignmentTransform( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, atomMap: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> object:
    """
    GetAlignmentTransform( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, atomMap: AtomPairsParameters = [], weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> object
        Compute the transformation required to align a molecule
             
              The 3D transformation required to align the specied conformation in the probe molecule
              to a specified conformation in the reference molecule is computed so that the root mean
              squared distance between a specified set of atoms is minimized
             
             ARGUMENTS
              - prbMol    molecule that is to be aligned
              - refMol    molecule used as the reference for the alignment
              - prbCid    ID of the conformation in the probe to be used 
                               for the alignment (defaults to first conformation)
              - refCid    ID of the conformation in the ref molecule to which 
                               the alignment is computed (defaults to first conformation)
              - atomMap   a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                               used to compute the alignments. If this mapping is 
                               not specified an attempt is made to generate one by
                               substructure matching
              - weights   Optionally specify weights for each of the atom pairs
              - reflect   if true reflect the conformation of the probe molecule
              - maxIters  maximum number of iterations used in minimizing the RMSD
               
              RETURNS
              a tuple of (RMSD value, transform matrix) 
            
        

        C++ signature :
            _object* GetAlignmentTransform(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]])
    """
def GetAllowNontetrahedralChirality() -> bool:
    """
    GetAllowNontetrahedralChirality() -> bool
        returns whether or not recognition of non-tetrahedral chirality from 3D structures is enabled

        C++ signature :
            bool GetAllowNontetrahedralChirality()
    """
def GetAngleDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float:
    """
    GetAngleDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float
        Returns the angle in degrees between atoms i, j, k
        

        C++ signature :
            double GetAngleDeg(RDKit::Conformer,unsigned int,unsigned int,unsigned int)
    """
def GetAngleRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float:
    """
    GetAngleRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float
        Returns the angle in radians between atoms i, j, k
        

        C++ signature :
            double GetAngleRad(RDKit::Conformer,unsigned int,unsigned int,unsigned int)
    """
def GetAtomAlias( atom: Atom) -> str:
    """
    GetAtomAlias( atom: Atom) -> str
        Returns the atom's MDL alias text

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAtomAlias(RDKit::Atom const*)
    """
def GetAtomFeatures( mol: Mol, atomid: int, addchiral: bool = False) -> list:
    """
    GetAtomFeatures( mol: Mol, atomid: int, addchiral: bool = False) -> list
        Returns the Atom Features vector

        C++ signature :
            boost::python::list GetAtomFeatures(RDKit::ROMol,int [,bool=False])
    """
def GetAtomMatch( featMatch: AtomPairsParameters, maxAts: int = 1024) -> object:
    """
    GetAtomMatch( featMatch: AtomPairsParameters, maxAts: int = 1024) -> object
        Returns an empty list if any of the features passed in share an atom.
         Otherwise a list of lists of atom indices is returned.
        

        C++ signature :
            boost::python::api::object GetAtomMatch(boost::python::api::object [,int=1024])
    """
def GetAtomPairAtomCode( atom: Atom, branchSubtract: int = 0, includeChirality: bool = False) -> int:
    """
    GetAtomPairAtomCode( atom: Atom, branchSubtract: int = 0, includeChirality: bool = False) -> int
        Returns the atom code (hash) for an atom

        C++ signature :
            unsigned int GetAtomPairAtomCode(RDKit::Atom const* [,unsigned int=0 [,bool=False]])
    """
def GetAtomPairAtomInvGen(includeChirality: bool = False) -> AtomInvariantsGenerator:
    """
    GetAtomPairAtomInvGen(includeChirality: bool = False) -> AtomInvariantsGenerator
        Get an atom pair atom-invariant generator
        
          ARGUMENTS:
            - includeChirality: if set, chirality will be taken into account for invariants
          RETURNS: AtomInvariantsGenerator
        
        

        C++ signature :
            RDKit::AtomInvariantsGenerator* GetAtomPairAtomInvGen([ bool=False])
    """
def GetAtomPairCode( atom1Code: int, atom2Code: int, distance: int, includeChirality: bool = False) -> int:
    """
    GetAtomPairCode( atom1Code: int, atom2Code: int, distance: int, includeChirality: bool = False) -> int
        Returns the atom-pair code (hash) for a pair of atoms separated by a certain number of bonds

        C++ signature :
            unsigned int GetAtomPairCode(unsigned int,unsigned int,unsigned int [,bool=False])
    """
def GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
    GetAtomPairFingerprint( mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect
        Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<int>* GetAtomPairFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]])
    """
def GetAtomPairGenerator(minDistance: int = 1, maxDistance: int = 30, includeChirality: bool = False, use2D: bool = True, countSimulation: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64:
    """
    GetAtomPairGenerator(minDistance: int = 1, maxDistance: int = 30, includeChirality: bool = False, use2D: bool = True, countSimulation: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64
        Get an atom pair fingerprint generator
        
          ARGUMENTS:
            - minDistance: minimum distance between atoms to be considered in a pair, default is 1 bond
            - maxDistance: maximum distance between atoms to be considered in a pair, default is maxPathLen-1 bonds
            - includeChirality: if set, chirality will be used in the atom  invariants, this is ignored if atomInvariantsGenerator is provided
            - use2D: if set, the 2D (topological) distance matrix  will be used
            - countSimulation:  if set, use count simulation while  generating the fingerprint
            - countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
            - fpSize: size of the generated fingerprint, does not affect the sparse versions
            - atomInvariantsGenerator: atom invariants to be used during fingerprint generation
        
        This generator supports the following AdditionalOutput types:
            - atomToBits: which bits each atom is involved in
            - atomCounts: how many bits each atom sets
            - bitInfoMap: map from bitId to (atomId, radius) pairs
        
          RETURNS: FingerprintGenerator
        
        

        C++ signature :
            RDKit::FingerprintGenerator<unsigned long>* GetAtomPairGenerator([ unsigned int=1 [,unsigned int=30 [,bool=False [,bool=True [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None]]]]]]]])
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
def GetBestAlignmentTransform( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, map: AtomPairsParameters = [], maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> object:
    """
    GetBestAlignmentTransform( prbMol: Mol, refMol: Mol, prbCid: int = -1, refCid: int = -1, map: AtomPairsParameters = [], maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = [], reflect: bool = False, maxIters: int = 50) -> object
        Compute the optimal RMS, transformation and atom map for aligning
              two molecules, taking symmetry into account. Molecule coordinates
              are left unaltered.
            
              This function will attempt to align all permutations of matching atom
              orders in both molecules, for some molecules it will lead to 'combinatorial
              explosion' especially if hydrogens are present.
              Use 'GetAlignmentTransform' to align molecules without changing the atom order.
            
             ARGUMENTS
              - prbMol      molecule that is to be aligned
              - refMol      molecule used as the reference for the alignment
              - prbCid      ID of the conformation in the probe to be used 
                            for the alignment (defaults to first conformation)
              - refCid      ID of the conformation in the ref molecule to which 
                            the alignment is computed (defaults to first conformation)
              - map:        (optional) a list of lists of (probeAtomId, refAtomId)
                            tuples with the atom-atom mappings of the two
                            molecules. If not provided, these will be generated
                            using a substructure search.
              - maxMatches  (optional) if atomMap is empty, this will be the max number of
                            matches found in a SubstructMatch().
              - symmetrizeConjugatedTerminalGroups (optional) if set, conjugated
                            terminal functional groups (like nitro or carboxylate)
                            will be considered symmetrically.
              - weights     Optionally specify weights for each of the atom pairs
              - reflect     if true reflect the conformation of the probe molecule
              - maxIters    maximum number of iterations used in minimizing the RMSD
               
              RETURNS
              a tuple of (RMSD value, best transform matrix, best atom map)
            
        

        C++ signature :
            _object* GetBestAlignmentTransform(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,boost::python::api::object=[] [,int=1000000 [,bool=True [,boost::python::api::object=[] [,bool=False [,unsigned int=50]]]]]]]])
    """
def GetBestRMS( prbMol: Mol, refMol: Mol, prbId: int = -1, refId: int = -1, map: AtomPairsParameters = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = []) -> float:
    """
    GetBestRMS( prbMol: Mol, refMol: Mol, prbId: int = -1, refId: int = -1, map: AtomPairsParameters = None, maxMatches: int = 1000000, symmetrizeConjugatedTerminalGroups: bool = True, weights: AtomPairsParameters = []) -> float
        Returns the optimal RMS for aligning two molecules, taking
               symmetry into account. As a side-effect, the probe molecule is
               left in the aligned state.
              
               Note:
               This function will attempt to align all permutations of matching atom
               orders in both molecules, for some molecules it will lead to
               'combinatorial explosion' especially if hydrogens are present.
               Use 'rdkit.Chem.AllChem.AlignMol' to align molecules without changing
               the atom order.
              
               ARGUMENTS
                - prbMol:      the molecule to be aligned to the reference
                - refMol:      the reference molecule
                - prbId:       (optional) probe conformation to use
                - refId:       (optional) reference conformation to use
                - map:         (optional) a list of lists of (probeAtomId,refAtomId)
                               tuples with the atom-atom mappings of the two
                               molecules. If not provided, these will be generated
                               using a substructure search.
                - maxMatches:  (optional) if map isn't specified, this will be
                               the max number of matches found in a SubstructMatch()
                - symmetrizeConjugatedTerminalGroups:  (optional) if set, conjugated
                               terminal functional groups (like nitro or carboxylate)
                               will be considered symmetrically
                - weights:     (optional) weights for mapping
               
              RETURNS
              The best RMSD found
            
        

        C++ signature :
            double GetBestRMS(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1 [,boost::python::api::object=None [,int=1000000 [,bool=True [,boost::python::api::object=[]]]]]]])
    """
def GetBondLength( conf: Conformer, iAtomId: int, jAtomId: int) -> float:
    """
    GetBondLength( conf: Conformer, iAtomId: int, jAtomId: int) -> float
        Returns the bond length in angstrom between atoms i, j
        

        C++ signature :
            double GetBondLength(RDKit::Conformer,unsigned int,unsigned int)
    """
def GetChemDrawRxnAdjustParams() -> AdjustQueryParameters:
    """
    GetChemDrawRxnAdjustParams() -> AdjustQueryParameters
        (deprecated, see MatchOnlyAtRgroupsAdjustParams)
            Returns the chemdraw style adjustment parameters for reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters GetChemDrawRxnAdjustParams()
    """
def GetConnectivityInvariants( mol: Mol, includeRingMembership: bool = True) -> list:
    """
    GetConnectivityInvariants( mol: Mol, includeRingMembership: bool = True) -> list
        Returns connectivity invariants (ECFP-like) for a molecule.

        C++ signature :
            boost::python::list GetConnectivityInvariants(RDKit::ROMol [,bool=True])
    """
def GetCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list:
    """
    GetCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
def GetCrippenO3A( prbMol: Mol, refMol: Mol, prbCrippenContribs: list = [], refCrippenContribs: list = [], prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> O3A:
    """
    GetCrippenO3A( prbMol: Mol, refMol: Mol, prbCrippenContribs: list = [], refCrippenContribs: list = [], prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> O3A
        Get an O3A object with atomMap and weights vectors to overlay
              the probe molecule onto the reference molecule based on
              Crippen logP atom contributions
             
             ARGUMENTS
              - prbMol                   molecule that is to be aligned
              - refMol                   molecule used as the reference for the alignment
              - prbCrippenContribs       Crippen atom contributions for the probe molecule
                                         as a list of (logp, mr) tuples, as returned
                                         by _CalcCrippenContribs()
              - refCrippenContribs       Crippen atom contributions for the reference molecule
                                         as a list of (logp, mr) tuples, as returned
                                         by _CalcCrippenContribs()
              - prbCid                   ID of the conformation in the probe to be used 
                                         for the alignment (defaults to first conformation)
              - refCid                   ID of the conformation in the ref molecule to which 
                                         the alignment is computed (defaults to first conformation)
              - reflect                  if true reflect the conformation of the probe molecule
                                         (defaults to false)
              - maxIters                 maximum number of iterations used in minimizing the RMSD
                                         (defaults to 50)
              - options                  least 2 significant bits encode accuracy
                                         (0: maximum, 3: minimum; defaults to 0)
                                         bit 3 triggers local optimization of the alignment
                                         (no computation of the cost matrix; defaults: off)
              - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                         which shall be used for the alignment (defaults to [])
              - constraintWeights        optionally specify weights for each of the constraints
                                         (weights default to 100.0)
               
              RETURNS
              The O3A object
            
        

        C++ signature :
            RDKit::MolAlign::PyO3A* GetCrippenO3A(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,boost::python::list=[] [,boost::python::list=[] [,int=-1 [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
def GetCrippenO3AForProbeConfs( prbMol: Mol, refMol: Mol, numThreads: int = 1, prbCrippenContribs: list = [], refCrippenContribs: list = [], refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> tuple:
    """
    GetCrippenO3AForProbeConfs( prbMol: Mol, refMol: Mol, numThreads: int = 1, prbCrippenContribs: list = [], refCrippenContribs: list = [], refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> tuple
        Get a vector of O3A objects for the overlay of all 
              the probe molecule's conformations onto the reference molecule based on
              MMFF atom types and charges
             
             ARGUMENTS
              - prbMol                   molecule that is to be aligned
              - refMol                   molecule used as the reference for the alignment
              - numThreads :             the number of threads to use, only has an effect if
                                         the RDKit was built with thread support (defaults to 1)
              - prbCrippenContribs       Crippen atom contributions for the probe molecule
                                         as a list of (logp, mr) tuples, as returned
                                         by _CalcCrippenContribs()
              - refCrippenContribs       Crippen atom contributions for the reference molecule
                                         as a list of (logp, mr) tuples, as returned
                                         by _CalcCrippenContribs()
              - refCid                   ID of the conformation in the ref molecule to which 
                                         the alignment is computed (defaults to first conformation)
              - reflect                  if true reflect the conformation of the probe molecule
                                         (defaults to false)
              - maxIters                 maximum number of iterations used in minimizing the RMSD
                                         (defaults to 50)
              - options                  least 2 significant bits encode accuracy
                                         (0: maximum, 3: minimum; defaults to 0)
                                         bit 3 triggers local optimization of the alignment
                                         (no computation of the cost matrix; defaults: off)
              - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                         which shall be used for the alignment (defaults to [])
              - constraintWeights        optionally specify weights for each of the constraints
                                         (weights default to 100.0)
               
              RETURNS
              A vector of O3A objects
            
        

        C++ signature :
            boost::python::tuple GetCrippenO3AForProbeConfs(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=1 [,boost::python::list=[] [,boost::python::list=[] [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
def GetDefaultAdjustParams() -> AdjustQueryParameters:
    """
    GetDefaultAdjustParams() -> AdjustQueryParameters
        Returns the default adjustment parameters for reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters GetDefaultAdjustParams()
    """
def GetDefaultPickleProperties() -> int:
    """
    GetDefaultPickleProperties() -> int
        Get the current global mol pickler options.

        C++ signature :
            unsigned int GetDefaultPickleProperties()
    """
def GetDihedralDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float:
    """
    GetDihedralDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float
        Returns the dihedral angle in degrees between atoms i, j, k, l
        

        C++ signature :
            double GetDihedralDeg(RDKit::Conformer,unsigned int,unsigned int,unsigned int,unsigned int)
    """
def GetDihedralRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float:
    """
    GetDihedralRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float
        Returns the dihedral angle in radians between atoms i, j, k, l
        

        C++ signature :
            double GetDihedralRad(RDKit::Conformer,unsigned int,unsigned int,unsigned int,unsigned int)
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
def GetErGFingerprint( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object:
    """
    GetErGFingerprint( mol: Mol, atomTypes: AtomPairsParameters = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> object
        Returns the ErG fingerprint vector for a molecule

        C++ signature :
            _object* GetErGFingerprint(RDKit::ROMol [,boost::python::api::object=0 [,double=0.3 [,int=1 [,int=15]]]])
    """
@typing.overload
def GetExperimentalTorsions( mol: Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> tuple:
    """
    GetExperimentalTorsions( mol: Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> tuple
        returns information about the bonds corresponding to experimental torsions

        C++ signature :
            boost::python::tuple GetExperimentalTorsions(RDKit::ROMol [,bool=True [,bool=False [,bool=True [,bool=True [,unsigned int=2 [,bool=False]]]]]])

        C++ signature :
            boost::python::tuple GetExperimentalTorsions(RDKit::ROMol,RDKit::DGeomHelpers::EmbedParameters)
    """
@typing.overload
def GetExperimentalTorsions( mol: Mol, embedParams: EmbedParameters) -> tuple:
    pass
def GetFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list:
    """
    GetFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
def GetFeatureInvariants( mol: Mol) -> list:
    """
    GetFeatureInvariants( mol: Mol) -> list
        Returns feature invariants (FCFP-like) for a molecule.

        C++ signature :
            boost::python::list GetFeatureInvariants(RDKit::ROMol)
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
def GetHashedAtomPairFingerprint( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
    GetHashedAtomPairFingerprint( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect
        Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<int>* GetHashedAtomPairFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False [,bool=True [,int=-1]]]]]]]]])
    """
def GetHashedAtomPairFingerprintAsBitVect( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> ExplicitBitVect:
    """
    GetHashedAtomPairFingerprintAsBitVect( mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> ExplicitBitVect
        Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect

        C++ signature :
            ExplicitBitVect* GetHashedAtomPairFingerprintAsBitVect(RDKit::ROMol [,unsigned int=2048 [,unsigned int=1 [,unsigned int=30 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,unsigned int=4 [,bool=False [,bool=True [,int=-1]]]]]]]]]])
    """
def GetHashedMorganFingerprint( mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect:
    """
    GetHashedMorganFingerprint( mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect
        Returns a hashed Morgan fingerprint for a molecule

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* GetHashedMorganFingerprint(RDKit::ROMol,unsigned int [,unsigned int=2048 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
def GetHashedTopologicalTorsionFingerprint( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
    GetHashedTopologicalTorsionFingerprint( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect
        Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<long>* GetHashedTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=2048 [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]]])
    """
def GetHashedTopologicalTorsionFingerprintAsBitVect( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False) -> ExplicitBitVect:
    """
    GetHashedTopologicalTorsionFingerprintAsBitVect( mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False) -> ExplicitBitVect
        Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect

        C++ signature :
            ExplicitBitVect* GetHashedTopologicalTorsionFingerprintAsBitVect(RDKit::ROMol [,unsigned int=2048 [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,unsigned int=4 [,bool=False]]]]]]])
    """
def GetMACCSKeysFingerprint( mol: Mol) -> ExplicitBitVect:
    """
    GetMACCSKeysFingerprint( mol: Mol) -> ExplicitBitVect
        Returns the MACCS keys for a molecule as an ExplicitBitVect

        C++ signature :
            ExplicitBitVect* GetMACCSKeysFingerprint(RDKit::ROMol)
    """
def GetMolFrags( mol: Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: AtomPairsParameters = None, fragsMolAtomMapping: AtomPairsParameters = None) -> tuple:
    """
    GetMolFrags( mol: Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: AtomPairsParameters = None, fragsMolAtomMapping: AtomPairsParameters = None) -> tuple
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
def GetMoleculeBoundsMatrix( mol: Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False) -> object:
    """
    GetMoleculeBoundsMatrix( mol: Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False) -> object
        Returns the distance bounds matrix for a molecule
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - set15bounds : set bounds for 1-5 atom distances based on 
                            topology (otherwise stop at 1-4s)
            - scaleVDW : scale down the sum of VDW radii when setting the 
                         lower bounds for atoms less than 5 bonds apart 
            - doTriangleSmoothing : run triangle smoothing on the bounds 
                         matrix before returning it 
         RETURNS:
        
            the bounds matrix as a Numeric array with lower bounds in 
            the lower triangle and upper bounds in the upper triangle
        
        

        C++ signature :
            _object* GetMoleculeBoundsMatrix(RDKit::ROMol {lvalue} [,bool=True [,bool=False [,bool=True [,bool=False]]]])
    """
def GetMorganAtomInvGen(includeRingMembership: bool = False) -> AtomInvariantsGenerator:
    """
    GetMorganAtomInvGen(includeRingMembership: bool = False) -> AtomInvariantsGenerator
        Get a morgan atom invariants generator
        
          ARGUMENTS:
            - includeRingMembership: if set, whether or not the atom is in a ring will be used in the invariant list
        
          RETURNS: AtomInvariantsGenerator
        
        

        C++ signature :
            RDKit::AtomInvariantsGenerator* GetMorganAtomInvGen([ bool=False])
    """
def GetMorganBondInvGen(useBondTypes: bool = True, useChirality: bool = False) -> BondInvariantsGenerator:
    """
    GetMorganBondInvGen(useBondTypes: bool = True, useChirality: bool = False) -> BondInvariantsGenerator
        Get a morgan bond invariants generator
        
          ARGUMENTS:
            - useBondTypes: if set, bond types will be included as a part of the bond invariants
            - useChirality: if set, chirality information will be included as a part of the bond invariants
        
          RETURNS: BondInvariantsGenerator
        
        

        C++ signature :
            RDKit::BondInvariantsGenerator* GetMorganBondInvGen([ bool=True [,bool=False]])
    """
def GetMorganFeatureAtomInvGen(patterns: AtomPairsParameters = None) -> AtomInvariantsGenerator:
    """
    GetMorganFeatureAtomInvGen(patterns: AtomPairsParameters = None) -> AtomInvariantsGenerator
        Get a morgan feature atom invariants generator
        
          ARGUMENTS:
            - patterns: if provided should contain the queries used to assign atom-types. if not provided, feature definitions adapted from reference: Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will be used for Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.
        
          RETURNS: AtomInvariantsGenerator
        
        

        C++ signature :
            RDKit::AtomInvariantsGenerator* GetMorganFeatureAtomInvGen([ boost::python::api::object {lvalue}=None])
    """
def GetMorganFingerprint( mol: Mol, radius: int, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, useCounts: bool = True, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect:
    """
    GetMorganFingerprint( mol: Mol, radius: int, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, useCounts: bool = True, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect
        Returns a Morgan fingerprint for a molecule

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* GetMorganFingerprint(RDKit::ROMol,unsigned int [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,bool=True [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
def GetMorganFingerprintAsBitVect( mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> ExplicitBitVect:
    """
    GetMorganFingerprintAsBitVect( mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> ExplicitBitVect
        Returns a Morgan fingerprint for a molecule as a bit vector

        C++ signature :
            ExplicitBitVect* GetMorganFingerprintAsBitVect(RDKit::ROMol,unsigned int [,unsigned int=2048 [,boost::python::api::object=[] [,boost::python::api::object=[] [,bool=False [,bool=True [,bool=False [,boost::python::api::object=None [,bool=False]]]]]]]])
    """
def GetMorganGenerator(radius: int = 3, countSimulation: bool = False, includeChirality: bool = False, useBondTypes: bool = True, onlyNonzeroInvariants: bool = False, includeRingMembership: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None, bondInvariantsGenerator: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> FingeprintGenerator64:
    """
    GetMorganGenerator(radius: int = 3, countSimulation: bool = False, includeChirality: bool = False, useBondTypes: bool = True, onlyNonzeroInvariants: bool = False, includeRingMembership: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None, bondInvariantsGenerator: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> FingeprintGenerator64
        Get a morgan fingerprint generator
        
          ARGUMENTS:
            - radius:  the number of iterations to grow the fingerprint
            - countSimulation: if set, use count simulation while generating the fingerprint
            - includeChirality: if set, chirality information will be added to the generated fingerprint
            - useBondTypes: if set, bond types will be included as a part of the default bond invariants
            - countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
            - fpSize: size of the generated fingerprint, does not affect the sparse versions
            - atomInvariantsGenerator: atom invariants to be used during fingerprint generation
        
        This generator supports the following AdditionalOutput types:
            - atomToBits: which bits each atom is the center of
            - atomCounts: how many bits each atom sets
            - bitInfoMap: map from bitId to (atomId1, radius) pairs
        
          RETURNS: FingerprintGenerator
        
        

        C++ signature :
            RDKit::FingerprintGenerator<unsigned long>* GetMorganGenerator([ unsigned int=3 [,bool=False [,bool=False [,bool=True [,bool=False [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None [,boost::python::api::object {lvalue}=None [,bool=False]]]]]]]]]]])
    """
def GetMostSubstitutedCoreMatch( mol: Mol, core: Mol, matches: AtomPairsParameters) -> object:
    """
    GetMostSubstitutedCoreMatch( mol: Mol, core: Mol, matches: AtomPairsParameters) -> object
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
def GetO3A( prbMol: Mol, refMol: Mol, prbPyMMFFMolProperties: AtomPairsParameters = None, refPyMMFFMolProperties: AtomPairsParameters = None, prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> O3A:
    """
    GetO3A( prbMol: Mol, refMol: Mol, prbPyMMFFMolProperties: AtomPairsParameters = None, refPyMMFFMolProperties: AtomPairsParameters = None, prbCid: int = -1, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> O3A
        Get an O3A object with atomMap and weights vectors to overlay
              the probe molecule onto the reference molecule based on
              MMFF atom types and charges
             
             ARGUMENTS
              - prbMol                   molecule that is to be aligned
              - refMol                   molecule used as the reference for the alignment
              - prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returned
                                         by SetupMMFFForceField()
              - refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returned
                                         by SetupMMFFForceField()
              - prbCid                   ID of the conformation in the probe to be used 
                                         for the alignment (defaults to first conformation)
              - refCid                   ID of the conformation in the ref molecule to which 
                                         the alignment is computed (defaults to first conformation)
              - reflect                  if true reflect the conformation of the probe molecule
                                         (defaults to false)
              - maxIters                 maximum number of iterations used in minimizing the RMSD
                                         (defaults to 50)
              - options                  least 2 significant bits encode accuracy
                                         (0: maximum, 3: minimum; defaults to 0)
                                         bit 3 triggers local optimization of the alignment
                                         (no computation of the cost matrix; defaults: off)
              - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                         which shall be used for the alignment (defaults to [])
              - constraintWeights        optionally specify weights for each of the constraints
                                         (weights default to 100.0)
               
              RETURNS
              The O3A object
            
        

        C++ signature :
            RDKit::MolAlign::PyO3A* GetO3A(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
def GetO3AForProbeConfs( prbMol: Mol, refMol: Mol, numThreads: int = 1, prbPyMMFFMolProperties: AtomPairsParameters = None, refPyMMFFMolProperties: AtomPairsParameters = None, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> tuple:
    """
    GetO3AForProbeConfs( prbMol: Mol, refMol: Mol, numThreads: int = 1, prbPyMMFFMolProperties: AtomPairsParameters = None, refPyMMFFMolProperties: AtomPairsParameters = None, refCid: int = -1, reflect: bool = False, maxIters: int = 50, options: int = 0, constraintMap: list = [], constraintWeights: list = []) -> tuple
        Get a vector of O3A objects for the overlay of all 
              the probe molecule's conformations onto the reference molecule based on
              MMFF atom types and charges
             
             ARGUMENTS
              - prbMol                   molecule that is to be aligned
              - refMol                   molecule used as the reference for the alignment
              - numThreads :             the number of threads to use, only has an effect if
                                         the RDKit was built with thread support (defaults to 1)
                                         If set to zero, the max supported by the system will be used.
              - prbPyMMFFMolProperties   PyMMFFMolProperties object for the probe molecule as returned
                                         by SetupMMFFForceField()
              - refPyMMFFMolProperties   PyMMFFMolProperties object for the reference molecule as returned
                                         by SetupMMFFForceField()
              - refCid                   ID of the conformation in the ref molecule to which 
                                         the alignment is computed (defaults to first conformation)
              - reflect                  if true reflect the conformation of the probe molecule
                                         (defaults to false)
              - maxIters                 maximum number of iterations used in minimizing the RMSD
                                         (defaults to 50)
              - options                  least 2 significant bits encode accuracy
                                         (0: maximum, 3: minimum; defaults to 0)
                                         bit 3 triggers local optimization of the alignment
                                         (no computation of the cost matrix; defaults: off)
              - constraintMap            a vector of pairs of atom IDs (probe AtomId, ref AtomId)
                                         which shall be used for the alignment (defaults to [])
              - constraintWeights        optionally specify weights for each of the constraints
                                         (weights default to 100.0)
               
              RETURNS
              A vector of O3A objects
            
        

        C++ signature :
            boost::python::tuple GetO3AForProbeConfs(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=1 [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,bool=False [,unsigned int=50 [,unsigned int=0 [,boost::python::list=[] [,boost::python::list=[]]]]]]]]]])
    """
def GetPeriodicTable() -> PeriodicTable:
    """
    GetPeriodicTable() -> PeriodicTable
        Returns the application's PeriodicTable instance.
        
        

        C++ signature :
            RDKit::PeriodicTable* GetPeriodicTable()
    """
def GetPreferCoordGen() -> bool:
    """
    GetPreferCoordGen() -> bool
        Return whether or not the CoordGen library is used for coordinate generation in the RDKit depiction library.

        C++ signature :
            bool GetPreferCoordGen()
    """
def GetRDKitAtomInvGen() -> AtomInvariantsGenerator:
    """
    GetRDKitAtomInvGen() -> AtomInvariantsGenerator
        Get an RDKit atom invariants generator
        
          RETURNS: AtomInvariantsGenerator
        
        

        C++ signature :
            RDKit::AtomInvariantsGenerator* GetRDKitAtomInvGen()
    """
def GetRDKitFPGenerator(minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, countSimulation: bool = False, countBounds: AtomPairsParameters = None, fpSize: int = 2048, numBitsPerFeature: int = 2, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64:
    """
    GetRDKitFPGenerator(minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, countSimulation: bool = False, countBounds: AtomPairsParameters = None, fpSize: int = 2048, numBitsPerFeature: int = 2, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64
        Get an RDKit fingerprint generator
        
          ARGUMENTS:
            - minPath: the minimum path length (in bonds) to be included
            - maxPath: the maximum path length (in bonds) to be included
            - useHs: toggles inclusion of Hs in paths (if the molecule has explicit Hs)
            - branchedPaths: toggles generation of branched subgraphs, not just linear paths
            - useBondOrder: toggles inclusion of bond orders in the path hashes
            - countSimulation:  if set, use count simulation while  generating the fingerprint
            - countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
            - fpSize: size of the generated fingerprint, does not affect the sparse versions
            - numBitsPerFeature: the number of bits set per path/subgraph found
            - atomInvariantsGenerator: atom invariants to be used during fingerprint generation
        
        This generator supports the following AdditionalOutput types:
            - atomToBits: which bits each atom is involved in
            - atomCounts: how many bits each atom sets
            - bitPaths: map from bitId to vectors of bond indices for the individual subgraphs
        
          RETURNS: FingerprintGenerator
        
        

        C++ signature :
            RDKit::FingerprintGenerator<unsigned long>* GetRDKitFPGenerator([ unsigned int=1 [,unsigned int=7 [,bool=True [,bool=True [,bool=True [,bool=False [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,unsigned int=2 [,boost::python::api::object {lvalue}=None]]]]]]]]]])
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
def GetSparseCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list:
    """
    GetSparseCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetSparseCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
def GetSparseFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list:
    """
    GetSparseFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetSparseFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
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
def GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
    GetTopologicalTorsionFingerprint( mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect
        Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect

        C++ signature :
            RDKit::SparseIntVect<long>* GetTopologicalTorsionFingerprint(RDKit::ROMol [,unsigned int=4 [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=0 [,bool=False]]]]])
    """
def GetTopologicalTorsionGenerator(includeChirality: bool = False, torsionAtomCount: int = 4, countSimulation: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64:
    """
    GetTopologicalTorsionGenerator(includeChirality: bool = False, torsionAtomCount: int = 4, countSimulation: bool = True, countBounds: AtomPairsParameters = None, fpSize: int = 2048, atomInvariantsGenerator: AtomPairsParameters = None) -> FingeprintGenerator64
        Get an atom pair fingerprint generator
        
          ARGUMENTS:
            - includeChirality: includeChirality argument for both the default atom invariants generator and the fingerprint arguments
            - torsionAtomCount: the number of atoms to include in the "torsions"
            - countSimulation:  if set, use count simulation while  generating the fingerprint
            - countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
            - fpSize: size of the generated fingerprint, does not affect the sparse versions
            - atomInvariantsGenerator: atom invariants to be used during fingerprint generation
        
        This generator supports the following AdditionalOutput types:
            - atomToBits: which bits each atom is involved in
            - atomCounts: how many bits each atom sets
            - bitPaths: map from bitId to vectors of atom indices
        
          RETURNS: FingerprintGenerator
        
        

        C++ signature :
            RDKit::FingerprintGenerator<unsigned long>* GetTopologicalTorsionGenerator([ bool=False [,unsigned int=4 [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None]]]]]])
    """
def GetUFFAngleBendParams( mol: Mol, idx1: int, idx2: int, idx3: int) -> object:
    """
    GetUFFAngleBendParams( mol: Mol, idx1: int, idx2: int, idx3: int) -> object
        Retrieves UFF angle bend parameters for atoms with indexes idx1, idx2, idx3 as a (ka, theta0) tuple, or None if no parameters could be found

        C++ signature :
            _object* GetUFFAngleBendParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int)
    """
def GetUFFBondStretchParams( mol: Mol, idx1: int, idx2: int) -> object:
    """
    GetUFFBondStretchParams( mol: Mol, idx1: int, idx2: int) -> object
        Retrieves UFF bond stretch parameters for atoms with indexes idx1, idx2 as a (kb, r0) tuple, or None if no parameters could be found

        C++ signature :
            _object* GetUFFBondStretchParams(RDKit::ROMol,unsigned int,unsigned int)
    """
def GetUFFInversionParams( mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
    """
    GetUFFInversionParams( mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object
        Retrieves UFF inversion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a K float value, or None if no parameters could be found

        C++ signature :
            _object* GetUFFInversionParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int,unsigned int)
    """
def GetUFFTorsionParams( mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
    """
    GetUFFTorsionParams( mol: Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object
        Retrieves UFF torsion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a V float value, or None if no parameters could be found

        C++ signature :
            _object* GetUFFTorsionParams(RDKit::ROMol,unsigned int,unsigned int,unsigned int,unsigned int)
    """
def GetUFFVdWParams( mol: Mol, idx1: int, idx2: int) -> object:
    """
    GetUFFVdWParams( mol: Mol, idx1: int, idx2: int) -> object
        Retrieves UFF van der Waals parameters for atoms with indexes idx1, idx2 as a (x_ij, D_ij) tuple, or None if no parameters could be found

        C++ signature :
            _object* GetUFFVdWParams(RDKit::ROMol,unsigned int,unsigned int)
    """
def GetUSR( mol: Mol, confId: int = -1) -> list:
    """
    GetUSR( mol: Mol, confId: int = -1) -> list
        Returns a USR descriptor for one conformer of a molecule

        C++ signature :
            boost::python::list GetUSR(RDKit::ROMol [,int=-1])
    """
def GetUSRCAT( mol: Mol, atomSelections: AtomPairsParameters = None, confId: int = -1) -> list:
    """
    GetUSRCAT( mol: Mol, atomSelections: AtomPairsParameters = None, confId: int = -1) -> list
        Returns a USRCAT descriptor for one conformer of a molecule

        C++ signature :
            boost::python::list GetUSRCAT(RDKit::ROMol [,boost::python::api::object=None [,int=-1]])
    """
def GetUSRDistributions( coords: AtomPairsParameters, points: AtomPairsParameters = None) -> list:
    """
    GetUSRDistributions( coords: AtomPairsParameters, points: AtomPairsParameters = None) -> list
        Returns the four USR distance distributions for a set of coordinates

        C++ signature :
            boost::python::list GetUSRDistributions(boost::python::api::object [,boost::python::api::object=None])
    """
def GetUSRDistributionsFromPoints( coords: AtomPairsParameters, points: AtomPairsParameters) -> list:
    """
    GetUSRDistributionsFromPoints( coords: AtomPairsParameters, points: AtomPairsParameters) -> list
        Returns the USR distance distributions for a set of coordinates and points

        C++ signature :
            boost::python::list GetUSRDistributionsFromPoints(boost::python::api::object,boost::python::api::object)
    """
def GetUSRFromDistributions( distances: AtomPairsParameters) -> list:
    """
    GetUSRFromDistributions( distances: AtomPairsParameters) -> list
        Returns the USR descriptor from a set of distance distributions

        C++ signature :
            boost::python::list GetUSRFromDistributions(boost::python::api::object)
    """
def GetUSRScore( descriptor1: AtomPairsParameters, descriptor2: AtomPairsParameters, weights: AtomPairsParameters = []) -> float:
    """
    GetUSRScore( descriptor1: AtomPairsParameters, descriptor2: AtomPairsParameters, weights: AtomPairsParameters = []) -> float
        Returns the USR score for two USR or USRCAT descriptors

        C++ signature :
            double GetUSRScore(boost::python::api::object,boost::python::api::object [,boost::python::api::object=[]])
    """
def GetUseLegacyStereoPerception() -> bool:
    """
    GetUseLegacyStereoPerception() -> bool
        returns whether or not the legacy stereo perception code is being used

        C++ signature :
            bool GetUseLegacyStereoPerception()
    """
def HCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* HCountEqualsQueryAtom(int [,bool=False])
    """
def HCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HCountGreaterQueryAtom(int [,bool=False])
    """
def HCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HCountLessQueryAtom(int [,bool=False])
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
def HasAgentTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    HasAgentTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool
        tests if the agents of a queryReaction are the same as those of a reaction

        C++ signature :
            bool HasAgentTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
def HasBitVectPropWithValueQueryAtom( propname: str, val: ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> QueryAtom:
    """
    HasBitVectPropWithValueQueryAtom( propname: str, val: ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified explicit bit vector value.  The Tolerance is the allowed Tanimoto difference

        C++ signature :
            RDKit::QueryAtom* HasBitVectPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,ExplicitBitVect [,bool=False [,float=0]])
    """
def HasBoolPropWithValueQueryAtom( propname: str, val: bool, negate: bool = False) -> QueryAtom:
    """
    HasBoolPropWithValueQueryAtom( propname: str, val: bool, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified boolean value.

        C++ signature :
            RDKit::QueryAtom* HasBoolPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
def HasBoolPropWithValueQueryBond( propname: str, val: bool, negate: bool = False) -> QueryBond:
    """
    HasBoolPropWithValueQueryBond( propname: str, val: bool, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified boolean value.

        C++ signature :
            RDKit::QueryBond* HasBoolPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
def HasChiralTagQueryAtom(negate: bool = False) -> QueryAtom:
    """
    HasChiralTagQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when HasChiralTag is True.

        C++ signature :
            RDKit::QueryAtom* HasChiralTagQueryAtom([ bool=False])
    """
def HasDoublePropWithValueQueryAtom( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryAtom:
    """
    HasDoublePropWithValueQueryAtom( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified value +- tolerance

        C++ signature :
            RDKit::QueryAtom* HasDoublePropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
def HasDoublePropWithValueQueryBond( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryBond:
    """
    HasDoublePropWithValueQueryBond( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified value +- tolerance

        C++ signature :
            RDKit::QueryBond* HasDoublePropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
def HasIntPropWithValueQueryAtom( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryAtom:
    """
    HasIntPropWithValueQueryAtom( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified int value.

        C++ signature :
            RDKit::QueryAtom* HasIntPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
def HasIntPropWithValueQueryBond( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryBond:
    """
    HasIntPropWithValueQueryBond( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified int value.

        C++ signature :
            RDKit::QueryBond* HasIntPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
def HasProductTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    HasProductTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool
        tests if the products of a queryReaction are substructures of the products of a reaction

        C++ signature :
            bool HasProductTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
def HasPropQueryAtom( propname: str, negate: bool = False) -> QueryAtom:
    """
    HasPropQueryAtom( propname: str, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' exists in the atom.

        C++ signature :
            RDKit::QueryAtom* HasPropQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasPropQueryBond( propname: str, negate: bool = False) -> QueryBond:
    """
    HasPropQueryBond( propname: str, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' exists in the bond.

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasReactantTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    HasReactantTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool
        tests if the reactants of a queryReaction are substructures of the reactants of a reaction

        C++ signature :
            bool HasReactantTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
def HasReactionAtomMapping( arg1: ChemicalReaction) -> bool:
    """
    HasReactionAtomMapping( arg1: ChemicalReaction) -> bool
        tests if a reaction obtains any atom mapping

        C++ signature :
            bool HasReactionAtomMapping(RDKit::ChemicalReaction)
    """
def HasReactionSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction, includeAgents: bool = False) -> bool:
    """
    HasReactionSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction, includeAgents: bool = False) -> bool
        tests if the queryReaction is a substructure of a reaction

        C++ signature :
            bool HasReactionSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction [,bool=False])
    """
def HasStringPropWithValueQueryAtom( propname: str, val: str, negate: bool = False) -> QueryAtom:
    """
    HasStringPropWithValueQueryAtom( propname: str, val: str, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified string value.

        C++ signature :
            RDKit::QueryAtom* HasStringPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasStringPropWithValueQueryBond( propname: str, val: str, negate: bool = False) -> QueryBond:
    """
    HasStringPropWithValueQueryBond( propname: str, val: str, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified string value.

        C++ signature :
            RDKit::QueryBond* HasStringPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HybridizationEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* HybridizationEqualsQueryAtom(int [,bool=False])
    """
def HybridizationGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HybridizationGreaterQueryAtom(int [,bool=False])
    """
def HybridizationLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HybridizationLessQueryAtom(int [,bool=False])
    """
def InNRingsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* InNRingsEqualsQueryAtom(int [,bool=False])
    """
def InNRingsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* InNRingsGreaterQueryAtom(int [,bool=False])
    """
def InNRingsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* InNRingsLessQueryAtom(int [,bool=False])
    """
def IsAliphaticQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsAliphaticQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsAliphatic is True.

        C++ signature :
            RDKit::QueryAtom* IsAliphaticQueryAtom([ bool=False])
    """
def IsAromaticQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsAromaticQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsAromatic is True.

        C++ signature :
            RDKit::QueryAtom* IsAromaticQueryAtom([ bool=False])
    """
def IsBridgeheadQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsBridgeheadQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsBridgehead is True.

        C++ signature :
            RDKit::QueryAtom* IsBridgeheadQueryAtom([ bool=False])
    """
def IsCoordGenSupportAvailable() -> bool:
    """
    IsCoordGenSupportAvailable() -> bool
        Returns whether RDKit was built with CoordGen support.

        C++ signature :
            bool IsCoordGenSupportAvailable()
    """
def IsInRingQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsInRingQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsInRing is True.

        C++ signature :
            RDKit::QueryAtom* IsInRingQueryAtom([ bool=False])
    """
def IsReactionTemplateMoleculeAgent( molecule: Mol, agentThreshold: float) -> bool:
    """
    IsReactionTemplateMoleculeAgent( molecule: Mol, agentThreshold: float) -> bool
        tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold

        C++ signature :
            bool IsReactionTemplateMoleculeAgent(RDKit::ROMol,double)
    """
def IsUnsaturatedQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsUnsaturatedQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsUnsaturated is True.

        C++ signature :
            RDKit::QueryAtom* IsUnsaturatedQueryAtom([ bool=False])
    """
def IsotopeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* IsotopeEqualsQueryAtom(int [,bool=False])
    """
def IsotopeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* IsotopeGreaterQueryAtom(int [,bool=False])
    """
def IsotopeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* IsotopeLessQueryAtom(int [,bool=False])
    """
def JSONToMols( jsonBlock: str, params: AtomPairsParameters = None) -> tuple:
    """
    JSONToMols( jsonBlock: str, params: AtomPairsParameters = None) -> tuple
        Convert JSON to a tuple of molecules
        
            ARGUMENTS:
              - jsonBlock: the molecule to work with
              - params: (optional) JSONParseParameters controlling the JSON parsing
            RETURNS:
              a tuple of Mols
        

        C++ signature :
            boost::python::tuple JSONToMols(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=None])
    """
def KDG() -> EmbedParameters:
    """
    KDG() -> EmbedParameters
        Returns an EmbedParameters object for the KDG method.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* KDG()
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
def LayeredFingerprint( mol: Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, branchedPaths: bool = True, fromAtoms: AtomPairsParameters = 0) -> ExplicitBitVect:
    """
    LayeredFingerprint( mol: Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, branchedPaths: bool = True, fromAtoms: AtomPairsParameters = 0) -> ExplicitBitVect
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
def LoadDefaultRingSystemTemplates() -> None:
    """
    LoadDefaultRingSystemTemplates() -> None
        Loads the default ring system templates and removes existing ones, if present.

        C++ signature :
            void LoadDefaultRingSystemTemplates()
    """
def MAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MAtom is True.

        C++ signature :
            RDKit::QueryAtom* MAtomQueryAtom([ bool=False])
    """
def MHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MHAtom is True.

        C++ signature :
            RDKit::QueryAtom* MHAtomQueryAtom([ bool=False])
    """
def MMFFGetMoleculeForceField( mol: Mol, pyMMFFMolProperties: MMFFMolProperties, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> ForceField:
    """
    MMFFGetMoleculeForceField( mol: Mol, pyMMFFMolProperties: MMFFMolProperties, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> ForceField
        returns a MMFF force field for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - pyMMFFMolProperties : PyMMFFMolProperties object as returned
                          by MMFFGetMoleculeProperties()
            - nonBondedThresh : used to exclude long-range non-bonded
                          interactions (defaults to 100.0)
            - confId : indicates which conformer to optimize
            - ignoreInterfragInteractions : if true, nonbonded terms between
                          fragments will not be added to the forcefield
        
        

        C++ signature :
            ForceFields::PyForceField* MMFFGetMoleculeForceField(RDKit::ROMol {lvalue},ForceFields::PyMMFFMolProperties* [,double=100.0 [,int=-1 [,bool=True]]])
    """
def MMFFGetMoleculeProperties( mol: Mol, mmffVariant: str = 'MMFF94', mmffVerbosity: int = 0) -> MMFFMolProperties:
    """
    MMFFGetMoleculeProperties( mol: Mol, mmffVariant: str = 'MMFF94', mmffVerbosity: int = 0) -> MMFFMolProperties
        returns a PyMMFFMolProperties object for a
          molecule, which is required by MMFFGetMoleculeForceField()
          and can be used to get/set MMFF properties
        
          
          ARGUMENTS:
        
            - mol : the molecule of interest
            - mmffVariant : "MMFF94" or "MMFF94s"
                          (defaults to "MMFF94")
            - mmffVerbosity : 0: none; 1: low; 2: high (defaults to 0).
        
        

        C++ signature :
            ForceFields::PyMMFFMolProperties* MMFFGetMoleculeProperties(RDKit::ROMol {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='MMFF94' [,unsigned int=0]])
    """
def MMFFHasAllMoleculeParams( mol: Mol) -> bool:
    """
    MMFFHasAllMoleculeParams( mol: Mol) -> bool
        checks if MMFF parameters are available for all of a molecule's atoms
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
        
        

        C++ signature :
            bool MMFFHasAllMoleculeParams(RDKit::ROMol)
    """
def MMFFOptimizeMolecule( mol: Mol, mmffVariant: str = 'MMFF94', maxIters: int = 200, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int:
    """
    MMFFOptimizeMolecule( mol: Mol, mmffVariant: str = 'MMFF94', maxIters: int = 200, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int
        uses MMFF to optimize a molecule's structure
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - mmffVariant : "MMFF94" or "MMFF94s"
            - maxIters : the maximum number of iterations (defaults to 200)
            - nonBondedThresh : used to exclude long-range non-bonded
                         interactions (defaults to 100.0)
            - confId : indicates which conformer to optimize
            - ignoreInterfragInteractions : if true, nonbonded terms between
                         fragments will not be added to the forcefield
        
         RETURNS: 0 if the optimization converged, -1 if the forcefield could
                  not be set up, 1 if more iterations are required.
        
        

        C++ signature :
            int MMFFOptimizeMolecule(RDKit::ROMol {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='MMFF94' [,int=200 [,double=100.0 [,int=-1 [,bool=True]]]]])
    """
def MMFFOptimizeMoleculeConfs( self: Mol, numThreads: int = 1, maxIters: int = 200, mmffVariant: str = 'MMFF94', nonBondedThresh: float = 100.0, ignoreInterfragInteractions: bool = True) -> object:
    """
    MMFFOptimizeMoleculeConfs( self: Mol, numThreads: int = 1, maxIters: int = 200, mmffVariant: str = 'MMFF94', nonBondedThresh: float = 100.0, ignoreInterfragInteractions: bool = True) -> object
        uses MMFF to optimize all of a molecule's conformations
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - numThreads : the number of threads to use, only has an effect if the RDKit
                           was built with thread support (defaults to 1)
                           If set to zero, the max supported by the system will be used.
            - maxIters : the maximum number of iterations (defaults to 200)
            - mmffVariant : "MMFF94" or "MMFF94s"
            - nonBondedThresh : used to exclude long-range non-bonded
                          interactions (defaults to 100.0)
            - ignoreInterfragInteractions : if true, nonbonded terms between
                          fragments will not be added to the forcefield.
        
        RETURNS: a list of (not_converged, energy) 2-tuples. 
            If not_converged is 0 the optimization converged for that conformer.
        
        

        C++ signature :
            boost::python::api::object MMFFOptimizeMoleculeConfs(RDKit::ROMol {lvalue} [,int=1 [,int=200 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='MMFF94' [,double=100.0 [,bool=True]]]]])
    """
def MMFFSanitizeMolecule( mol: Mol) -> int:
    """
    MMFFSanitizeMolecule( mol: Mol) -> int
        sanitizes a molecule according to MMFF requirements.
        
            - mol : the molecule of interest.
        
        

        C++ signature :
            unsigned int MMFFSanitizeMolecule(RDKit::ROMol {lvalue})
    """
def MQNs_( mol: Mol, force: bool = False) -> list:
    """
    MQNs_( mol: Mol, force: bool = False) -> list

        C++ signature :
            boost::python::list MQNs_(RDKit::ROMol [,bool=False])
    """
def MakePropertyRangeQuery( name: str, min: float, max: float) -> PropertyRangeQuery:
    """
    MakePropertyRangeQuery( name: str, min: float, max: float) -> PropertyRangeQuery
        Generates a Range property for the specified property, between min and max
        query = MakePropertyRangeQuery('exactmw', 0, 500)
        query.Match( mol )

        C++ signature :
            Queries::RangeQuery<double, RDKit::ROMol const&, true>* MakePropertyRangeQuery(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double,double)
    """
def MassEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* MassEqualsQueryAtom(int [,bool=False])
    """
def MassGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MassGreaterQueryAtom(int [,bool=False])
    """
def MassLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MassLessQueryAtom(int [,bool=False])
    """
def MatchOnlyAtRgroupsAdjustParams() -> AdjustQueryParameters:
    """
    MatchOnlyAtRgroupsAdjustParams() -> AdjustQueryParameters
        Only match at the specified rgroup locations in the reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters MatchOnlyAtRgroupsAdjustParams()
    """
def MergeQueryHs( mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> Mol:
    """
    MergeQueryHs( mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> Mol
        merges hydrogens into their neighboring atoms as queries

        C++ signature :
            RDKit::ROMol* MergeQueryHs(RDKit::ROMol [,bool=False [,bool=False]])
    """
def MetadataFromPNGFile( filename: AtomPairsParameters) -> dict:
    """
    MetadataFromPNGFile( filename: AtomPairsParameters) -> dict
        Returns a dict with all metadata from the PNG file. Keys are strings, values are bytes.

        C++ signature :
            boost::python::dict MetadataFromPNGFile(boost::python::api::object)
    """
def MetadataFromPNGString( png: AtomPairsParameters) -> dict:
    """
    MetadataFromPNGString( png: AtomPairsParameters) -> dict
        Returns a dict with all metadata from the PNG string. Keys are strings, values are bytes.

        C++ signature :
            boost::python::dict MetadataFromPNGString(boost::python::api::object)
    """
def MinRingSizeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* MinRingSizeEqualsQueryAtom(int [,bool=False])
    """
def MinRingSizeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MinRingSizeGreaterQueryAtom(int [,bool=False])
    """
def MinRingSizeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MinRingSizeLessQueryAtom(int [,bool=False])
    """
def MissingChiralTagQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MissingChiralTagQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MissingChiralTag is True.

        C++ signature :
            RDKit::QueryAtom* MissingChiralTagQueryAtom([ bool=False])
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
def MolFragmentToCXSmarts( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, isomericSmarts: bool = True) -> str:
    """
    MolFragmentToCXSmarts( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, isomericSmarts: bool = True) -> str
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
def MolFragmentToCXSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0) -> str:
    """
    MolFragmentToCXSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0) -> str
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
def MolFragmentToCXSmiles( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    pass
def MolFragmentToSmarts( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, isomericSmarts: bool = True) -> str:
    """
    MolFragmentToSmarts( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, isomericSmarts: bool = True) -> str
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
def MolFragmentToSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0) -> str:
    """
    MolFragmentToSmiles( mol: Mol, params: SmilesWriteParams, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0) -> str
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
def MolFragmentToSmiles( mol: Mol, atomsToUse: AtomPairsParameters, bondsToUse: AtomPairsParameters = 0, atomSymbols: AtomPairsParameters = 0, bondSymbols: AtomPairsParameters = 0, isomericSmiles: bool = True, kekuleSmiles: bool = False, rootedAtAtom: int = -1, canonical: bool = True, allBondsExplicit: bool = False, allHsExplicit: bool = False) -> str:
    pass
def MolFromFASTA( text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0) -> Mol:
    """
    MolFromFASTA( text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0) -> Mol
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
def MolFromHELM( text: AtomPairsParameters, sanitize: bool = True) -> Mol:
    """
    MolFromHELM( text: AtomPairsParameters, sanitize: bool = True) -> Mol
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
def MolFromMolBlock( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol:
    """
    MolFromMolBlock( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> Mol
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
def MolFromMrvBlock( mrvBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> Mol:
    """
    MolFromMrvBlock( mrvBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> Mol
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
def MolFromPDBBlock( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol:
    """
    MolFromPDBBlock( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True, flavor: int = 0, proximityBonding: bool = True) -> Mol
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
def MolFromPNGFile( filename: str, params: AtomPairsParameters = None) -> Mol:
    """
    MolFromPNGFile( filename: str, params: AtomPairsParameters = None) -> Mol
        Construct a molecule from metadata in a PNG file.
        
             ARGUMENTS:
        
               - filename: the PNG filename
        
               - params: used to provide optional parameters for the metadata parsing
        
             RETURNS:
               a Mol object, None on failure.

        C++ signature :
            RDKit::ROMol* MolFromPNGFile(char const* [,boost::python::api::object=None])
    """
def MolFromPNGString( png: AtomPairsParameters, params: AtomPairsParameters = None) -> Mol:
    """
    MolFromPNGString( png: AtomPairsParameters, params: AtomPairsParameters = None) -> Mol
        Construct a molecule from metadata in a PNG string.
        
             ARGUMENTS:
        
               - png: the PNG string
        
               - params: used to provide optional parameters for the metadata parsing
        
             RETURNS:
               a Mol object, None on failure.
          

        C++ signature :
            RDKit::ROMol* MolFromPNGString(boost::python::api::object [,boost::python::api::object=None])
    """
def MolFromQuerySLN( SLN: str, mergeHs: bool = True, debugParser: bool = False) -> Mol:
    """
    MolFromQuerySLN( SLN: str, mergeHs: bool = True, debugParser: bool = False) -> Mol
        Construct a query molecule from an SLN string.
        
          ARGUMENTS:
        
            - SLN: the SLN string
        
            - mergeHs: (optional) toggles the merging of explicit Hs in the query into the attached
              heavy atoms. Defaults to False.
        
          RETURNS:
        
            a Mol object suitable for using in substructure queries, None on failure.
        
        

        C++ signature :
            RDKit::ROMol* MolFromQuerySLN(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
def MolFromRDKitSVG( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> Mol:
    """
    MolFromRDKitSVG( molBlock: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> Mol
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
def MolFromSLN( SLN: str, sanitize: bool = True, debugParser: bool = False) -> Mol:
    """
    MolFromSLN( SLN: str, sanitize: bool = True, debugParser: bool = False) -> Mol
        Construct a molecule from an SLN string.
        
            ARGUMENTS:
        
            - SLN: the SLN string
        
            - sanitize: (optional) toggles sanitization of the molecule.
              Defaults to True.
        
          RETURNS:
        
            a Mol object, None on failure.
        
          NOTE: the SLN should not contain query information or properties. To build a
            query from SLN, use MolFromQuerySLN.
        
        

        C++ signature :
            RDKit::ROMol* MolFromSLN(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
def MolFromSequence( text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0) -> Mol:
    """
    MolFromSequence( text: AtomPairsParameters, sanitize: bool = True, flavor: int = 0) -> Mol
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
def MolFromSmarts( SMARTS: AtomPairsParameters, mergeHs: bool = False, replacements: dict = {}) -> Mol:
    """
    MolFromSmarts( SMARTS: AtomPairsParameters, mergeHs: bool = False, replacements: dict = {}) -> Mol
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
def MolFromSmarts( SMARTS: AtomPairsParameters, params: SmartsParserParams) -> Mol:
    pass
@typing.overload
def MolFromSmiles( SMILES: AtomPairsParameters, params: SmilesParserParams) -> Mol:
    """
    MolFromSmiles( SMILES: AtomPairsParameters, params: SmilesParserParams) -> Mol
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
def MolFromSmiles( SMILES: AtomPairsParameters, sanitize: bool = True, replacements: dict = {}) -> Mol:
    pass
def MolFromTPLBlock( tplBlock: AtomPairsParameters, sanitize: bool = True, skipFirstConf: bool = False) -> Mol:
    """
    MolFromTPLBlock( tplBlock: AtomPairsParameters, sanitize: bool = True, skipFirstConf: bool = False) -> Mol
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
def MolFromXYZBlock( xyzFileName: AtomPairsParameters) -> Mol:
    """
    MolFromXYZBlock( xyzFileName: AtomPairsParameters) -> Mol
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
def MolMetadataToPNGFile( mol: Mol, filename: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object:
    """
    MolMetadataToPNGFile( mol: Mol, filename: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object
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
def MolMetadataToPNGString( mol: Mol, png: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object:
    """
    MolMetadataToPNGString( mol: Mol, png: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeMol: bool = False) -> object
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
def MolToCXSmiles( mol: Mol, params: SmilesWriteParams, flags: int = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL) -> str:
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
def MolToJSON( mol: Mol, params: AtomPairsParameters = None) -> str:
    """
    MolToJSON( mol: Mol, params: AtomPairsParameters = None) -> str
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
def MolsFromCDXML( cdxml: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> tuple:
    """
    MolsFromCDXML( cdxml: AtomPairsParameters, sanitize: bool = True, removeHs: bool = True) -> tuple
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
def MolsFromPNGFile( filename: str, tag: str = 'rdkitPKL', params: AtomPairsParameters = None) -> object:
    """
    MolsFromPNGFile( filename: str, tag: str = 'rdkitPKL', params: AtomPairsParameters = None) -> object
        returns a tuple of molecules constructed from the PNG file

        C++ signature :
            boost::python::api::object MolsFromPNGFile(char const* [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='rdkitPKL' [,boost::python::api::object=None]])
    """
def MolsFromPNGString( png: AtomPairsParameters, tag: str = 'rdkitPKL', params: AtomPairsParameters = None) -> tuple:
    """
    MolsFromPNGString( png: AtomPairsParameters, tag: str = 'rdkitPKL', params: AtomPairsParameters = None) -> tuple
        returns a tuple of molecules constructed from the PNG string

        C++ signature :
            boost::python::tuple MolsFromPNGString(boost::python::api::object [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='rdkitPKL' [,boost::python::api::object=None]])
    """
def MolsToJSON( mols: AtomPairsParameters, params: AtomPairsParameters = None) -> str:
    """
    MolsToJSON( mols: AtomPairsParameters, params: AtomPairsParameters = None) -> str
        Convert a set of molecules to JSON
        
            ARGUMENTS:
              - mols: the molecules to work with
            RETURNS:
              a string
        

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolsToJSON(boost::python::api::object [,boost::python::api::object=None])
    """
def MrvBlockIsReaction( mrvData: str) -> bool:
    """
    MrvBlockIsReaction( mrvData: str) -> bool
        returns whether or not an MRV block contains reaction data

        C++ signature :
            bool MrvBlockIsReaction(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def MrvFileIsReaction( filename: str) -> bool:
    """
    MrvFileIsReaction( filename: str) -> bool
        returns whether or not an MRV file contains reaction data

        C++ signature :
            bool MrvFileIsReaction(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def MurckoDecompose( mol: Mol) -> Mol:
    """
    MurckoDecompose( mol: Mol) -> Mol
        Do a Murcko decomposition and return the scaffold

        C++ signature :
            RDKit::ROMol* MurckoDecompose(RDKit::ROMol)
    """
def NonHydrogenDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeEqualsQueryAtom(int [,bool=False])
    """
def NonHydrogenDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeGreaterQueryAtom(int [,bool=False])
    """
def NonHydrogenDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeLessQueryAtom(int [,bool=False])
    """
def NormalizeDepiction( mol: Mol, confId: int = -1, canonicalize: int = 1, scaleFactor: float = -1.0) -> float:
    """
    NormalizeDepiction( mol: Mol, confId: int = -1, canonicalize: int = 1, scaleFactor: float = -1.0) -> float
        Normalizes the 2D depiction.
        If canonicalize is != 0, the depiction is subjected to a canonical
        transformation such that its main axis is aligned along the X axis
        (canonicalize >0, the default) or the Y axis (canonicalize <0).
        If canonicalize is 0, no canonicalization takes place.
        If scaleFactor is <0.0 (the default) the depiction is scaled such
        that bond lengths conform to RDKit standards. The applied scaling
        factor is returned.
        
        ARGUMENTS:
        
        mol          - the molecule to be normalized
        confId       - (optional) the id of the reference conformation to use
        canonicalize - (optional) if != 0, a canonical transformation is
                       applied: if >0 (the default), the main molecule axis is
                       aligned to the X axis, if <0 to the Y axis.
                       If 0, no canonical transformation is applied.
        scaleFactor  - (optional) if >0.0, the scaling factor to apply. The default
                       (-1.0) means that the depiction is automatically scaled
                       such that bond lengths are the standard RDKit ones.
        
        RETURNS: the applied scaling factor.

        C++ signature :
            double NormalizeDepiction(RDKit::ROMol {lvalue} [,int=-1 [,int=1 [,double=-1.0]]])
    """
def NumAliphaticHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
def NumAliphaticHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
def NumAliphaticHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsEqualsQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsGreaterQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsLessQueryAtom(int [,bool=False])
    """
def OptimizeMolecule( ff: ForceField, maxIters: int = 200) -> int:
    """
    OptimizeMolecule( ff: ForceField, maxIters: int = 200) -> int
        uses the supplied force field to optimize a molecule's structure
        
         
         ARGUMENTS:
        
            - ff : the force field
            - maxIters : the maximum number of iterations (defaults to 200)
        
         RETURNS: 0 if the optimization converged, 1 if more iterations are required.
        
        

        C++ signature :
            int OptimizeMolecule(ForceFields::PyForceField {lvalue} [,int=200])
    """
def OptimizeMoleculeConfs( mol: Mol, ff: ForceField, numThreads: int = 1, maxIters: int = 200) -> object:
    """
    OptimizeMoleculeConfs( mol: Mol, ff: ForceField, numThreads: int = 1, maxIters: int = 200) -> object
        uses the supplied force field to optimize all of a molecule's conformations
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - ff : the force field
            - numThreads : the number of threads to use, only has an effect if the RDKit
                           was built with thread support (defaults to 1)
                           If set to zero, the max supported by the system will be used.
            - maxIters : the maximum number of iterations (defaults to 200)
        
         RETURNS: a list of (not_converged, energy) 2-tuples. 
             If not_converged is 0 the optimization converged for that conformer.
        
        

        C++ signature :
            boost::python::api::object OptimizeMoleculeConfs(RDKit::ROMol {lvalue},ForceFields::PyForceField {lvalue} [,int=1 [,int=200]])
    """
def PEOE_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    PEOE_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list PEOE_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
def ParseMolQueryDefFile( fileobj: AtomPairsParameters, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict:
    """
    ParseMolQueryDefFile( fileobj: AtomPairsParameters, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict
        reads query definitions from a simply formatted file
        

        C++ signature :
            boost::python::dict ParseMolQueryDefFile(boost::python::api::object {lvalue} [,bool=True [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='\t' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='//' [,unsigned int=0 [,unsigned int=1]]]]])
    """
def PathToSubmol( mol: Mol, path: AtomPairsParameters, useQuery: bool = False, atomMap: AtomPairsParameters = None) -> Mol:
    """
    PathToSubmol( mol: Mol, path: AtomPairsParameters, useQuery: bool = False, atomMap: AtomPairsParameters = None) -> Mol

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
def PreprocessReaction( reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue') -> object:
    """
    PreprocessReaction( reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue') -> object
        A function for preprocessing reactions with more specific queries.
        Queries are indicated by labels on atoms (molFileAlias property by default)
        When these labels are found, more specific queries are placed on the atoms.
        By default, the available quieries come from 
          FilterCatalog.GetFlattenedFunctionalGroupHierarchy(True)n
        Sample Usage:
          >>> from rdkit import Chem, RDConfig
          >>> from rdkit.Chem import MolFromSmiles, AllChem
          >>> from rdkit.Chem.rdChemReactions import PreprocessReaction
          >>> import os
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','boronic1.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          >>> nWarn
          0
          >>> nError
          0
          >>> nReacts
          2
          >>> nProds
          1
          >>> reactantLabels
          (((0, 'halogen.bromine.aromatic'),), ((1, 'boronicacid'),))
        
        If there are functional group labels in the input reaction (via atoms with molFileValue properties),
        the corresponding atoms will have queries added to them so that they only match such things. We can
        see this here:
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> r1 = rxn.GetReactantTemplate(0)
          >>> m1 = Chem.MolFromSmiles('CCBr')
          >>> m2 = Chem.MolFromSmiles('c1ccccc1Br')
          
        These both match because the reaction file itself just has R1-Br:
          >>> m1.HasSubstructMatch(r1)
          True
          >>> m2.HasSubstructMatch(r1)
          True
        
        After preprocessing, we only match the aromatic Br:
          >>> d = PreprocessReaction(rxn)
          >>> m1.HasSubstructMatch(r1)
          False
          >>> m2.HasSubstructMatch(r1)
          True
        
        We also support or queries in the values field (separated by commas):
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','azide_reaction.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> reactantLabels = PreprocessReaction(rxn)[-1]
          >>> reactantLabels
          (((1, 'azide'),), ((1, 'carboxylicacid,acidchloride'),))
          >>> m1 = Chem.MolFromSmiles('CC(=O)O')
          >>> m2 = Chem.MolFromSmiles('CC(=O)Cl')
          >>> m3 = Chem.MolFromSmiles('CC(=O)N')
          >>> r2 = rxn.GetReactantTemplate(1)
          >>> m1.HasSubstructMatch(r2)
          True
          >>> m2.HasSubstructMatch(r2)
          True
          >>> m3.HasSubstructMatch(r2)
          False
        
        unrecognized final group types are returned as None:
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value1.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'boromicacid'
        
        One unrecognized group type in a comma-separated list makes the whole thing fail:
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value2.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'carboxylicacid,acidchlroide'
          >>> testFile = os.path.join(RDConfig.RDCodeDir,'Chem','SimpleEnum','test_data','bad_value3.rxn')
          >>> rxn = AllChem.ReactionFromRxnFile(testFile)
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          Traceback (most recent call last):
            ...
          KeyError: 'carboxyliccaid,acidchloride'
          >>> rxn = rdChemReactions.ChemicalReaction()
          >>> rxn.Initialize()
          >>> nWarn,nError,nReacts,nProds,reactantLabels = PreprocessReaction(rxn)
          >>> reactantLabels
          ()
          >>> reactantLabels == ()
          True
        

        C++ signature :
            boost::python::api::object PreprocessReaction(RDKit::ChemicalReaction {lvalue} [,boost::python::dict={} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='molFileValue']])
    """
def QAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    QAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when QAtom is True.

        C++ signature :
            RDKit::QueryAtom* QAtomQueryAtom([ bool=False])
    """
def QHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    QHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when QHAtom is True.

        C++ signature :
            RDKit::QueryAtom* QHAtomQueryAtom([ bool=False])
    """
def RDKFingerprint( mol: Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: AtomPairsParameters = 0, fromAtoms: AtomPairsParameters = 0, atomBits: AtomPairsParameters = None, bitInfo: AtomPairsParameters = None) -> ExplicitBitVect:
    """
    RDKFingerprint( mol: Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: AtomPairsParameters = 0, fromAtoms: AtomPairsParameters = 0, atomBits: AtomPairsParameters = None, bitInfo: AtomPairsParameters = None) -> ExplicitBitVect
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
def RandomTransform( mol: Mol, cid: int = -1, seed: int = -1) -> None:
    """
    RandomTransform( mol: Mol, cid: int = -1, seed: int = -1) -> None
        Perform a random transformation on a molecule
             
             ARGUMENTS
              - mol    molecule that is to be transformed
              - cid    ID of the conformation in the mol to be transformed
                       (defaults to first conformation)
              - seed   seed used to initialize the random generator
                       (defaults to -1, that is no seeding)
               
            
        

        C++ signature :
            void RandomTransform(RDKit::ROMol {lvalue} [,int=-1 [,int=-1]])
    """
def ReactionFromMolecule( arg1: Mol) -> ChemicalReaction:
    """
    ReactionFromMolecule( arg1: Mol) -> ChemicalReaction
        construct a ChemicalReaction from an molecule if the RXN role property of the molecule is set

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromMolecule(RDKit::ROMol)
    """
def ReactionFromMrvBlock( rxnblock: AtomPairsParameters, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction:
    """
    ReactionFromMrvBlock( rxnblock: AtomPairsParameters, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction
        construct a ChemicalReaction from a string in Marvin (mrv) format

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromMrvBlock(boost::python::api::object [,bool=False [,bool=False]])

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromMrvBlock(boost::python::api::object [,bool=False [,bool=False]])
    """
def ReactionFromMrvFile( filename: str, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction:
    """
    ReactionFromMrvFile( filename: str, sanitize: bool = False, removeHs: bool = False) -> ChemicalReaction
        construct a ChemicalReaction from an Marvin (mrv) rxn file

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromMrvFile(char const* [,bool=False [,bool=False]])
    """
def ReactionFromPNGFile( arg1: str) -> ChemicalReaction:
    """
    ReactionFromPNGFile( arg1: str) -> ChemicalReaction
        construct a ChemicalReaction from metadata in a PNG file

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromPNGFile(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def ReactionFromPNGString( arg1: str) -> ChemicalReaction:
    """
    ReactionFromPNGString( arg1: str) -> ChemicalReaction
        construct a ChemicalReaction from an string with PNG data

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromPNGString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def ReactionFromRxnBlock( rxnblock: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction:
    """
    ReactionFromRxnBlock( rxnblock: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction
        construct a ChemicalReaction from a string in MDL rxn format

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromRxnBlock(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False [,bool=True]]])
    """
def ReactionFromRxnFile( filename: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction:
    """
    ReactionFromRxnFile( filename: str, sanitize: bool = False, removeHs: bool = False, strictParsing: bool = True) -> ChemicalReaction
        construct a ChemicalReaction from an MDL rxn file

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromRxnFile(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False [,bool=True]]])
    """
def ReactionFromSmarts( SMARTS: str, replacements: dict = {}, useSmiles: bool = False) -> ChemicalReaction:
    """
    ReactionFromSmarts( SMARTS: str, replacements: dict = {}, useSmiles: bool = False) -> ChemicalReaction
        construct a ChemicalReaction from a reaction SMARTS string. 
        see the documentation for rdkit.Chem.MolFromSmiles for an explanation
        of the replacements argument.

        C++ signature :
            RDKit::ChemicalReaction* ReactionFromSmarts(char const* [,boost::python::dict={} [,bool=False]])
    """
def ReactionMetadataToPNGFile( mol: ChemicalReaction, filename: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeMol: bool = False) -> object:
    """
    ReactionMetadataToPNGFile( mol: ChemicalReaction, filename: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeMol: bool = False) -> object
        Reads the contents of a PNG file and adds metadata about a reaction to it. The modified file contents are returned.

        C++ signature :
            boost::python::api::object ReactionMetadataToPNGFile(RDKit::ChemicalReaction,boost::python::api::object [,bool=True [,bool=True [,bool=False [,bool=False]]]])
    """
def ReactionMetadataToPNGString( mol: ChemicalReaction, pngdata: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeRxn: bool = False) -> object:
    """
    ReactionMetadataToPNGString( mol: ChemicalReaction, pngdata: AtomPairsParameters, includePkl: bool = True, includeSmiles: bool = True, includeSmarts: bool = False, includeRxn: bool = False) -> object
        Adds metadata about a reaction to the PNG string passed in.The modified string is returned.

        C++ signature :
            boost::python::api::object ReactionMetadataToPNGString(RDKit::ChemicalReaction,boost::python::api::object [,bool=True [,bool=True [,bool=False [,bool=False]]]])
    """
def ReactionToMolecule( reaction: ChemicalReaction) -> Mol:
    """
    ReactionToMolecule( reaction: ChemicalReaction) -> Mol
        construct a molecule for a ChemicalReaction with RXN role property set

        C++ signature :
            RDKit::ROMol* ReactionToMolecule(RDKit::ChemicalReaction)
    """
def ReactionToMrvBlock( reaction: ChemicalReaction, prettyPrint: bool = False) -> str:
    """
    ReactionToMrvBlock( reaction: ChemicalReaction, prettyPrint: bool = False) -> str
        construct a string in Marvin (MRV) rxn format for a ChemicalReaction

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToMrvBlock(RDKit::ChemicalReaction [,bool=False])
    """
def ReactionToMrvFile( reaction: ChemicalReaction, filename: str, prettyPrint: bool = False) -> None:
    """
    ReactionToMrvFile( reaction: ChemicalReaction, filename: str, prettyPrint: bool = False) -> None
        write a Marvin (MRV) rxn file for a ChemicalReaction

        C++ signature :
            void ReactionToMrvFile(RDKit::ChemicalReaction,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def ReactionToRxnBlock( reaction: ChemicalReaction, separateAgents: bool = False, forceV3000: bool = False) -> str:
    """
    ReactionToRxnBlock( reaction: ChemicalReaction, separateAgents: bool = False, forceV3000: bool = False) -> str
        construct a string in MDL rxn format for a ChemicalReaction

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToRxnBlock(RDKit::ChemicalReaction [,bool=False [,bool=False]])
    """
def ReactionToSmarts( reaction: ChemicalReaction) -> str:
    """
    ReactionToSmarts( reaction: ChemicalReaction) -> str
        construct a reaction SMARTS string for a ChemicalReaction

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToSmarts(RDKit::ChemicalReaction)
    """
def ReactionToSmiles( reaction: ChemicalReaction, canonical: bool = True) -> str:
    """
    ReactionToSmiles( reaction: ChemicalReaction, canonical: bool = True) -> str
        construct a reaction SMILES string for a ChemicalReaction

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToSmiles(RDKit::ChemicalReaction [,bool=True])
    """
def ReactionToV3KRxnBlock( reaction: ChemicalReaction, separateAgents: bool = False) -> str:
    """
    ReactionToV3KRxnBlock( reaction: ChemicalReaction, separateAgents: bool = False) -> str
        construct a string in MDL v3000 rxn format for a ChemicalReaction

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > ReactionToV3KRxnBlock(RDKit::ChemicalReaction [,bool=False])
    """
def ReactionsFromCDXMLBlock( rxnblock: AtomPairsParameters, sanitize: bool = False, removeHs: bool = False) -> object:
    """
    ReactionsFromCDXMLBlock( rxnblock: AtomPairsParameters, sanitize: bool = False, removeHs: bool = False) -> object
        construct a tuple of ChemicalReactions from a string in CDXML format

        C++ signature :
            boost::python::api::object ReactionsFromCDXMLBlock(boost::python::api::object [,bool=False [,bool=False]])
    """
def ReactionsFromCDXMLFile( filename: str, sanitize: bool = False, removeHs: bool = False) -> object:
    """
    ReactionsFromCDXMLFile( filename: str, sanitize: bool = False, removeHs: bool = False) -> object
        construct a tuple of ChemicalReactions from a CDXML rxn file

        C++ signature :
            boost::python::api::object ReactionsFromCDXMLFile(char const* [,bool=False [,bool=False]])
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
def ReduceProductToSideChains( product: Mol, addDummyAtoms: bool = True) -> Mol:
    """
    ReduceProductToSideChains( product: Mol, addDummyAtoms: bool = True) -> Mol
        reduce the product of a reaction to the side chains added by the reaction.              The output is a molecule with attached wildcards indicating where the product was attached.              The dummy atom has the same reaction-map number as the product atom (if available).

        C++ signature :
            RDKit::ROMol* ReduceProductToSideChains(boost::shared_ptr<RDKit::ROMol> [,bool=True])
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
def RemoveMappingNumbersFromReactions( reaction: ChemicalReaction) -> None:
    """
    RemoveMappingNumbersFromReactions( reaction: ChemicalReaction) -> None
        Removes the mapping numbers from the molecules of a reaction

        C++ signature :
            void RemoveMappingNumbersFromReactions(RDKit::ChemicalReaction)
    """
def RemoveStereochemistry( mol: Mol) -> None:
    """
    RemoveStereochemistry( mol: Mol) -> None
        Removes all stereochemistry info from the molecule.
        
        

        C++ signature :
            void RemoveStereochemistry(RDKit::ROMol {lvalue})
    """
def RenumberAtoms( mol: Mol, newOrder: AtomPairsParameters) -> Mol:
    """
    RenumberAtoms( mol: Mol, newOrder: AtomPairsParameters) -> Mol
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
def ReplaceCore( mol: Mol, core: Mol, matches: AtomPairsParameters, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> Mol:
    """
    ReplaceCore( mol: Mol, core: Mol, matches: AtomPairsParameters, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> Mol
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
def RingBondCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* RingBondCountEqualsQueryAtom(int [,bool=False])
    """
def RingBondCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* RingBondCountGreaterQueryAtom(int [,bool=False])
    """
def RingBondCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* RingBondCountLessQueryAtom(int [,bool=False])
    """
def SMR_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    SMR_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list SMR_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
def SanitizeMol( mol: Mol, sanitizeOps: int = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors: bool = False) -> SanitizeFlags:
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
def SanitizeRxn( rxn: ChemicalReaction, sanitizeOps: int = 4294967295, params: AdjustQueryParameters = AdjustQueryParameters(), catchErrors: bool = False) -> SanitizeFlags:
    """
    SanitizeRxn( rxn: ChemicalReaction, sanitizeOps: int = 4294967295, params: AdjustQueryParameters = AdjustQueryParameters(), catchErrors: bool = False) -> SanitizeFlags
        Does some sanitization of the reactant and product templates of a reaction.
        
            - The reaction is modified in place.
            - If sanitization fails, an exception will be thrown unless catchErrors is set
        
          ARGUMENTS:
        
            - rxn: the reaction to be modified
            - sanitizeOps: (optional) reaction sanitization operations to be carried out
              these should be constructed by or'ing together the
              operations in rdkit.Chem.rdChemReactions.SanitizeFlags
            - optional adjustment parameters for changing the meaning of the substructure
              matching done in the templates.  The default is 
              rdkit.Chem.rdChemReactions.DefaultRxnAdjustParams which aromatizes
              kekule structures if possible.
            - catchErrors: (optional) if provided, instead of raising an exception
              when sanitization fails (the default behavior), the 
              first operation that failed (as defined in rdkit.Chem.rdChemReactions.SanitizeFlags)
              is returned. Zero is returned on success.
        
          The operations carried out by default are:
            1) fixRGroups(): sets R group labels on mapped dummy atoms when possible
            2) fixAtomMaps(): attempts to set atom maps on unmapped R groups
            3) adjustTemplate(): calls adjustQueryProperties() on all reactant templates
            4) fixHs(): merges explicit Hs in the reactant templates that don't map to heavy atoms
        

        C++ signature :
            RDKit::RxnOps::SanitizeRxnFlags SanitizeRxn(RDKit::ChemicalReaction {lvalue} [,unsigned long=4294967295 [,RDKit::MolOps::AdjustQueryParameters=<rdkit.Chem.rdmolops.AdjustQueryParameters object at 0x7f2c64d90360> [,bool=False]]])
    """
def SetAllowNontetrahedralChirality( arg1: bool) -> None:
    """
    SetAllowNontetrahedralChirality( arg1: bool) -> None
        toggles recognition of non-tetrahedral chirality from 3D structures

        C++ signature :
            void SetAllowNontetrahedralChirality(bool)
    """
def SetAngleDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None:
    """
    SetAngleDeg( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None
        Sets the angle in degrees between atoms i, j, k; all atoms bonded to atom k are moved
        

        C++ signature :
            void SetAngleDeg(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,double)
    """
def SetAngleRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None:
    """
    SetAngleRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None
        Sets the angle in radians between atoms i, j, k; all atoms bonded to atom k are moved
        

        C++ signature :
            void SetAngleRad(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,double)
    """
def SetAromaticity( mol: Mol, model: AromaticityModel = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT) -> None:
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
def SetBondLength( conf: Conformer, iAtomId: int, jAtomId: int, value: float) -> None:
    """
    SetBondLength( conf: Conformer, iAtomId: int, jAtomId: int, value: float) -> None
        Sets the bond length in angstrom between atoms i, j; all atoms bonded to atom j are moved
        

        C++ signature :
            void SetBondLength(RDKit::Conformer {lvalue},unsigned int,unsigned int,double)
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
def SetDihedralDeg( arg1: Conformer, arg2: int, arg3: int, arg4: int, arg5: int, arg6: float) -> None:
    """
    SetDihedralDeg( arg1: Conformer, arg2: int, arg3: int, arg4: int, arg5: int, arg6: float) -> None
        Sets the dihedral angle in degrees between atoms i, j, k, l; all atoms bonded to atom l are moved
        

        C++ signature :
            void SetDihedralDeg(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,unsigned int,double)
    """
def SetDihedralRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int, value: float) -> None:
    """
    SetDihedralRad( conf: Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int, value: float) -> None
        Sets the dihedral angle in radians between atoms i, j, k, l; all atoms bonded to atom l are moved
        

        C++ signature :
            void SetDihedralRad(RDKit::Conformer {lvalue},unsigned int,unsigned int,unsigned int,unsigned int,double)
    """
def SetDoubleBondNeighborDirections( mol: Mol, conf: AtomPairsParameters = None) -> None:
    """
    SetDoubleBondNeighborDirections( mol: Mol, conf: AtomPairsParameters = None) -> None
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
def SetPreferCoordGen( val: bool) -> None:
    """
    SetPreferCoordGen( val: bool) -> None
        Sets whether or not the CoordGen library should be preferred to the RDKit depiction library.

        C++ signature :
            void SetPreferCoordGen(bool)
    """
def SetRingSystemTemplates( templatePath: str) -> None:
    """
    SetRingSystemTemplates( templatePath: str) -> None
        Loads the ring system templates from the specified file to be used in 2D coordinate generation. Each template must be a single line in the file represented using CXSMILES, and the structure should be a single ring system. Throws a DepictException if any templates are invalid.

        C++ signature :
            void SetRingSystemTemplates(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
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
def ShapeProtrudeDist( mol1: Mol, mol2: Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True, allowReordering: bool = True) -> float:
    """
    ShapeProtrudeDist( mol1: Mol, mol2: Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True, allowReordering: bool = True) -> float
        Compute the shape protrude distance between two molecule based on a predefined alignment
          
          ARGUMENTS:
            - mol1 : The first molecule of interest 
            - mol2 : The second molecule of interest 
            - confId1 : Conformer in the first molecule (defaults to first conformer) 
            - confId2 : Conformer in the second molecule (defaults to first conformer) 
            - gridSpacing : resolution of the grid used to encode the molecular shapes 
            - bitsPerPoint : number of bit used to encode the occupancy at each grid point 
                                  defaults to two bits per grid point 
            - vdwScale : Scaling factor for the radius of the atoms to determine the base radius 
                        used in the encoding - grid points inside this sphere carry the maximum occupancy 
            - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased 
                         from layer to layer from the maximum value 
            - maxLayers : the maximum number of layers - defaults to the number of bits 
                          used per grid point - e.g. two bits per grid point will allow 3 layers 
            - ignoreHs : when set, the contribution of Hs to the shape will be ignored
            - allowReordering : when set, the order will be automatically updated so that the value calculated
                                is the protrusion of the smaller shape from the larger one.
        

        C++ signature :
            double ShapeProtrudeDist(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,double=0.5 [,RDKit::DiscreteValueVect::DiscreteValueType=rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE [,double=0.8 [,double=0.25 [,int=-1 [,bool=True [,bool=True]]]]]]]]])
    """
def ShapeTanimotoDist( mol1: Mol, mol2: Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float:
    """
    ShapeTanimotoDist( mol1: Mol, mol2: Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float
        Compute the shape tanimoto distance between two molecule based on a predefined alignment
          
          ARGUMENTS:
            - mol1 : The first molecule of interest 
            - mol2 : The second molecule of interest 
            - confId1 : Conformer in the first molecule (defaults to first conformer) 
            - confId2 : Conformer in the second molecule (defaults to first conformer) 
            - gridSpacing : resolution of the grid used to encode the molecular shapes 
            - bitsPerPoint : number of bits used to encode the occupancy at each grid point 
                                  defaults to two bits per grid point 
            - vdwScale : Scaling factor for the radius of the atoms to determine the base radius 
                        used in the encoding - grid points inside this sphere carry the maximum occupancy 
            - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased 
                         from layer to layer from the maximum value 
            - maxLayers : the maximum number of layers - defaults to the number of bits 
                          used per grid point - e.g. two bits per grid point will allow 3 layers 
            - ignoreHs : when set, the contribution of Hs to the shape will be ignored
        

        C++ signature :
            double ShapeTanimotoDist(RDKit::ROMol,RDKit::ROMol [,int=-1 [,int=-1 [,double=0.5 [,RDKit::DiscreteValueVect::DiscreteValueType=rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE [,double=0.8 [,double=0.25 [,int=-1 [,bool=True]]]]]]]])
    """
def ShapeTverskyIndex( mol1: Mol, mol2: Mol, alpha: float, beta: float, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float:
    """
    ShapeTverskyIndex( mol1: Mol, mol2: Mol, alpha: float, beta: float, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float
        Compute the shape tversky index between two molecule based on a predefined alignment
          
          ARGUMENTS:
            - mol1 : The first molecule of interest 
            - mol2 : The second molecule of interest 
            - alpha : first parameter of the Tversky index
            - beta : second parameter of the Tversky index
            - confId1 : Conformer in the first molecule (defaults to first conformer) 
            - confId2 : Conformer in the second molecule (defaults to first conformer) 
            - gridSpacing : resolution of the grid used to encode the molecular shapes 
            - bitsPerPoint : number of bits used to encode the occupancy at each grid point 
                                  defaults to two bits per grid point 
            - vdwScale : Scaling factor for the radius of the atoms to determine the base radius 
                        used in the encoding - grid points inside this sphere carry the maximum occupancy 
            - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased 
                         from layer to layer from the maximum value 
            - maxLayers : the maximum number of layers - defaults to the number of bits 
                          used per grid point - e.g. two bits per grid point will allow 3 layers 
            - ignoreHs : when set, the contribution of Hs to the shape will be ignored
        

        C++ signature :
            double ShapeTverskyIndex(RDKit::ROMol,RDKit::ROMol,double,double [,int=-1 [,int=-1 [,double=0.5 [,RDKit::DiscreteValueVect::DiscreteValueType=rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE [,double=0.8 [,double=0.25 [,int=-1 [,bool=True]]]]]]]])
    """
def SlogP_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    SlogP_VSA_( mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list

        C++ signature :
            boost::python::list SlogP_VSA_(RDKit::ROMol [,boost::python::api::object=[] [,bool=False]])
    """
def SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier:
    """
    SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier

        C++ signature :
            RDKit::SmilesMolSupplier* SmilesMolSupplierFromText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
    """
def SortMatchesByDegreeOfCoreSubstitution( mol: Mol, core: Mol, matches: AtomPairsParameters) -> object:
    """
    SortMatchesByDegreeOfCoreSubstitution( mol: Mol, core: Mol, matches: AtomPairsParameters) -> object
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
def SplitMolByPDBChainId( mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False) -> dict:
    """
    SplitMolByPDBChainId( mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False) -> dict
        Splits a molecule into pieces based on PDB chain information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by chain id with molecules as the values
        
        

        C++ signature :
            boost::python::dict SplitMolByPDBChainId(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
def SplitMolByPDBResidues( mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False) -> dict:
    """
    SplitMolByPDBResidues( mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False) -> dict
        Splits a molecule into pieces based on PDB residue information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by residue name with molecules as the values
        
        

        C++ signature :
            boost::python::dict SplitMolByPDBResidues(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
def StraightenDepiction( mol: Mol, confId: int = -1, minimizeRotation: bool = False) -> None:
    """
    StraightenDepiction( mol: Mol, confId: int = -1, minimizeRotation: bool = False) -> None
        Rotate the 2D depiction such that the majority of bonds have a
          30-degree angle with the X axis.
          ARGUMENTS:
        
          mol              - the molecule to be rotated.
          confId           - (optional) the id of the reference conformation to use.
          minimizeRotation - (optional) if False (the default), the molecule
                             is rotated such that the majority of bonds have an angle
                             with the X axis of 30 or 90 degrees. If True, the minimum
                             rotation is applied such that the majority of bonds have
                             an angle with the X axis of 0, 30, 60, or 90 degrees,
                             with the goal of altering the initial orientation as
                             little as possible .

        C++ signature :
            void StraightenDepiction(RDKit::ROMol {lvalue} [,int=-1 [,bool=False]])
    """
def TotalDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* TotalDegreeEqualsQueryAtom(int [,bool=False])
    """
def TotalDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalDegreeGreaterQueryAtom(int [,bool=False])
    """
def TotalDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalDegreeLessQueryAtom(int [,bool=False])
    """
def TotalValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* TotalValenceEqualsQueryAtom(int [,bool=False])
    """
def TotalValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalValenceGreaterQueryAtom(int [,bool=False])
    """
def TotalValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalValenceLessQueryAtom(int [,bool=False])
    """
def TransformConformer( arg1: Conformer, arg2: AtomPairsParameters) -> None:
    """
    TransformConformer( arg1: Conformer, arg2: AtomPairsParameters) -> None
        Transform the coordinates of a conformer

        C++ signature :
            void TransformConformer(RDKit::Conformer {lvalue},boost::python::api::object)
    """
def TranslateChiralFlagToStereoGroups( mol: Mol, zeroFlagGroupType: StereoGroupType = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND) -> None:
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
def UFFGetMoleculeForceField( mol: Mol, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> ForceField:
    """
    UFFGetMoleculeForceField( mol: Mol, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> ForceField
        returns a UFF force field for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - vdwThresh : used to exclude long-range van der Waals interactions
                          (defaults to 10.0)
            - confId : indicates which conformer to optimize
            - ignoreInterfragInteractions : if true, nonbonded terms between
                          fragments will not be added to the forcefield.
        
        

        C++ signature :
            ForceFields::PyForceField* UFFGetMoleculeForceField(RDKit::ROMol {lvalue} [,double=10.0 [,int=-1 [,bool=True]]])
    """
def UFFHasAllMoleculeParams( mol: Mol) -> bool:
    """
    UFFHasAllMoleculeParams( mol: Mol) -> bool
        checks if UFF parameters are available for all of a molecule's atoms
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest.
        
        

        C++ signature :
            bool UFFHasAllMoleculeParams(RDKit::ROMol)
    """
def UFFOptimizeMolecule( self: Mol, maxIters: int = 200, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int:
    """
    UFFOptimizeMolecule( self: Mol, maxIters: int = 200, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int
        uses UFF to optimize a molecule's structure
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - maxIters : the maximum number of iterations (defaults to 200)
            - vdwThresh : used to exclude long-range van der Waals interactions
                          (defaults to 10.0)
            - confId : indicates which conformer to optimize
            - ignoreInterfragInteractions : if true, nonbonded terms between
                          fragments will not be added to the forcefield.
        
         RETURNS: 0 if the optimization converged, 1 if more iterations are required.
        
        

        C++ signature :
            int UFFOptimizeMolecule(RDKit::ROMol {lvalue} [,int=200 [,double=10.0 [,int=-1 [,bool=True]]]])
    """
def UFFOptimizeMoleculeConfs( self: Mol, numThreads: int = 1, maxIters: int = 200, vdwThresh: float = 10.0, ignoreInterfragInteractions: bool = True) -> object:
    """
    UFFOptimizeMoleculeConfs( self: Mol, numThreads: int = 1, maxIters: int = 200, vdwThresh: float = 10.0, ignoreInterfragInteractions: bool = True) -> object
        uses UFF to optimize all of a molecule's conformations
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - numThreads : the number of threads to use, only has an effect if the RDKit
                           was built with thread support (defaults to 1)
                           If set to zero, the max supported by the system will be used.
            - maxIters : the maximum number of iterations (defaults to 200)
            - vdwThresh : used to exclude long-range van der Waals interactions
                          (defaults to 10.0)
            - ignoreInterfragInteractions : if true, nonbonded terms between
                          fragments will not be added to the forcefield.
        
         RETURNS: a list of (not_converged, energy) 2-tuples. 
             If not_converged is 0 the optimization converged for that conformer.
        
        

        C++ signature :
            boost::python::api::object UFFOptimizeMoleculeConfs(RDKit::ROMol {lvalue} [,int=1 [,int=200 [,double=10.0 [,bool=True]]]])
    """
def UnfoldedRDKFingerprintCountBased( mol: Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: AtomPairsParameters = 0, fromAtoms: AtomPairsParameters = 0, atomBits: AtomPairsParameters = None, bitInfo: AtomPairsParameters = None) -> ULongSparseIntVect:
    """
    UnfoldedRDKFingerprintCountBased( mol: Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: AtomPairsParameters = 0, fromAtoms: AtomPairsParameters = 0, atomBits: AtomPairsParameters = None, bitInfo: AtomPairsParameters = None) -> ULongSparseIntVect
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
def UpdateProductsStereochemistry( reaction: ChemicalReaction) -> None:
    """
    UpdateProductsStereochemistry( reaction: ChemicalReaction) -> None
        Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.          This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.          The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction

        C++ signature :
            void UpdateProductsStereochemistry(RDKit::ChemicalReaction*)
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
def XAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    XAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when XAtom is True.

        C++ signature :
            RDKit::QueryAtom* XAtomQueryAtom([ bool=False])
    """
def XHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    XHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when XHAtom is True.

        C++ signature :
            RDKit::QueryAtom* XHAtomQueryAtom([ bool=False])
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
def molzipFragments( mols: AtomPairsParameters, params: MolzipParams = MolzipParams()) -> Mol:
    """
    molzipFragments( mols: AtomPairsParameters, params: MolzipParams = MolzipParams()) -> Mol
        zip together multiple molecules from an R group decomposition 
        using the given matching parameters.  The first molecule in the list
        must be the core

        C++ signature :
            RDKit::ROMol* molzipFragments(boost::python::api::object {lvalue} [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7f2c653ca0f0>])
    """
def srETKDGv3() -> EmbedParameters:
    """
    srETKDGv3() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 3 (small rings).

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* srETKDGv3()
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
AtomPairFP = rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP
AtomProps = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
BAD_DOUBLE_BOND_STEREO = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
BondProps = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
CHECK_CHIRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
CHECK_TETRAHEDRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
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
ETK_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
FINAL_CENTER_IN_VOLUME = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
FINAL_CHIRAL_BOUNDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
FIRST_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
INCHI_AVAILABLE = True
INITIAL_COORDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
KEKULE_ALL = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
LINEAR_DOUBLE_BOND = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
LayeredFingerprint_substructLayers = 7
MINIMIZE_FOURTH_DIMENSION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
MolProps = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
MorganFP = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP
NoConformers = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
RDKitFP = rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP
SANITIZE_ADJUSTHS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
SANITIZE_ADJUST_REACTANTS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
SANITIZE_ALL = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
SANITIZE_ATOM_MAPS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
SANITIZE_CLEANUP = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
SANITIZE_CLEANUPCHIRALITY = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
SANITIZE_CLEANUP_ORGANOMETALLICS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
SANITIZE_FINDRADICALS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
SANITIZE_KEKULIZE = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
SANITIZE_MERGEHS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
SANITIZE_NONE = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
SANITIZE_PROPERTIES = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
SANITIZE_RGROUP_NAMES = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
SANITIZE_SETAROMATICITY = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
SANITIZE_SETCONJUGATION = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
SANITIZE_SETHYBRIDIZATION = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
SANITIZE_SYMMRINGS = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
STEREO_ABSOLUTE = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
TopologicalTorsionFP = rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP
UNCONSTRAINED_ANIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
logger: rdkit.RDLogger.logger
templDir = '/scratch/toscopa1/src/rdkit/Data/'
