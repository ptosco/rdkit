from __future__ import annotations
import rdkit.Chem.rdFingerprintGenerator
import typing
import Boost.Python

__all__ = [
    "AdditionalOutput",
    "AtomInvariantsGenerator",
    "AtomPairFP",
    "AtomPairFingerprintOptions",
    "BondInvariantsGenerator",
    "FPType",
    "FingeprintGenerator32",
    "FingeprintGenerator64",
    "FingerprintOptions",
    "GetAtomPairAtomInvGen",
    "GetAtomPairGenerator",
    "GetCountFPs",
    "GetFPs",
    "GetMorganAtomInvGen",
    "GetMorganBondInvGen",
    "GetMorganFeatureAtomInvGen",
    "GetMorganGenerator",
    "GetRDKitAtomInvGen",
    "GetRDKitFPGenerator",
    "GetSparseCountFPs",
    "GetSparseFPs",
    "GetTopologicalTorsionGenerator",
    "MorganFP",
    "MorganFingerprintOptions",
    "RDKitFP",
    "RDKitFingerprintOptions",
    "TopologicalTorsionFP",
    "TopologicalTorsionFingerprintOptions"
]


class AdditionalOutput(Boost.Python.instance):
    @staticmethod
    def AllocateAtomCounts( arg1: AdditionalOutput) -> None: 
        """
        AllocateAtomCounts( arg1: AdditionalOutput) -> None
            synonym for CollectAtomCounts()

            C++ signature :
                void AllocateAtomCounts(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def AllocateAtomToBits( arg1: AdditionalOutput) -> None: 
        """
        AllocateAtomToBits( arg1: AdditionalOutput) -> None
            synonym for CollectAtomToBits()

            C++ signature :
                void AllocateAtomToBits(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def AllocateBitInfoMap( arg1: AdditionalOutput) -> None: 
        """
        AllocateBitInfoMap( arg1: AdditionalOutput) -> None
            synonym for CollectBitInfoMap()

            C++ signature :
                void AllocateBitInfoMap(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def AllocateBitPaths( arg1: AdditionalOutput) -> None: 
        """
        AllocateBitPaths( arg1: AdditionalOutput) -> None
            synonym for CollectBitPaths()

            C++ signature :
                void AllocateBitPaths(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def CollectAtomCounts( arg1: AdditionalOutput) -> None: 
        """
        CollectAtomCounts( arg1: AdditionalOutput) -> None
            toggles collection of information about the number of bits each atom is involved in

            C++ signature :
                void CollectAtomCounts(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def CollectAtomToBits( arg1: AdditionalOutput) -> None: 
        """
        CollectAtomToBits( arg1: AdditionalOutput) -> None
            toggle collection of information mapping each atom to the bits it is involved in.

            C++ signature :
                void CollectAtomToBits(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def CollectBitInfoMap( arg1: AdditionalOutput) -> None: 
        """
        CollectBitInfoMap( arg1: AdditionalOutput) -> None
            toggles collection of information mapping each atom to more detail about the atom environment (not available from all fingerprints)

            C++ signature :
                void CollectBitInfoMap(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def CollectBitPaths( arg1: AdditionalOutput) -> None: 
        """
        CollectBitPaths( arg1: AdditionalOutput) -> None
            toggles collection of information matching each atom to information about the paths it is involved in (not available from all fingerprints).

            C++ signature :
                void CollectBitPaths(RDKit::AdditionalOutput {lvalue})
        """
    @staticmethod
    def GetAtomCounts( arg1: AdditionalOutput) -> object: 
        """
        GetAtomCounts( arg1: AdditionalOutput) -> object

            C++ signature :
                boost::python::api::object GetAtomCounts(RDKit::AdditionalOutput)
        """
    @staticmethod
    def GetAtomToBits( arg1: AdditionalOutput) -> object: 
        """
        GetAtomToBits( arg1: AdditionalOutput) -> object

            C++ signature :
                boost::python::api::object GetAtomToBits(RDKit::AdditionalOutput)
        """
    @staticmethod
    def GetBitInfoMap( arg1: AdditionalOutput) -> object: 
        """
        GetBitInfoMap( arg1: AdditionalOutput) -> object

            C++ signature :
                boost::python::api::object GetBitInfoMap(RDKit::AdditionalOutput)
        """
    @staticmethod
    def GetBitPaths( arg1: AdditionalOutput) -> object: 
        """
        GetBitPaths( arg1: AdditionalOutput) -> object

            C++ signature :
                boost::python::api::object GetBitPaths(RDKit::AdditionalOutput)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 88
    pass
class AtomInvariantsGenerator(Boost.Python.instance):
    pass
class FingerprintOptions(Boost.Python.instance):
    @staticmethod
    def SetCountBounds( arg1: FingerprintOptions, arg2: AtomPairsParameters) -> None: 
        """
        SetCountBounds( arg1: FingerprintOptions, arg2: AtomPairsParameters) -> None
            set the bins for the count bounds

            C++ signature :
                void SetCountBounds(RDKit::FingerprintArguments {lvalue},boost::python::api::object)
        """
    @property
    def countSimulation(self) -> None:
        """
        use count simulation

        :type: None
        """
    @property
    def fpSize(self) -> None:
        """
        size of the fingerprints created

        :type: None
        """
    @property
    def includeChirality(self) -> None:
        """
        include chirality in atom invariants (not for all fingerprints)

        :type: None
        """
    @property
    def numBitsPerFeature(self) -> None:
        """
        number of bits to set for each feature

        :type: None
        """
    pass
class BondInvariantsGenerator(Boost.Python.instance):
    pass
class FPType(Boost.Python.enum, int):
    AtomPairFP = rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP
    MorganFP = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP
    RDKitFP = rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP
    TopologicalTorsionFP = rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP
    __slots__ = ()
    names = {'RDKitFP': rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP, 'MorganFP': rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP, 'AtomPairFP': rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP, 'TopologicalTorsionFP': rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP}
    values = {2: rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP, 1: rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP, 0: rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP, 3: rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP}
    pass
class FingeprintGenerator32(Boost.Python.instance):
    @staticmethod
    def GetCountFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect: 
        """
        GetCountFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect
            Generates a count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseIntVect containing fingerprint
            
            

            C++ signature :
                RDKit::SparseIntVect<unsigned int>* GetCountFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetCountFingerprintAsNumPy( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object: 
        """
        GetCountFingerprintAsNumPy( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object
            Generates a count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a numpy array containing the fingerprint
            
            

            C++ signature :
                boost::python::api::object GetCountFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ExplicitBitVect: 
        """
        GetFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ExplicitBitVect
            Generates a fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a ExplicitBitVect containing fingerprint
            
            

            C++ signature :
                ExplicitBitVect* GetFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetFingerprintAsNumPy( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object: 
        """
        GetFingerprintAsNumPy( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object
            Generates a fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a numpy array containing the fingerprint
            
            

            C++ signature :
                boost::python::api::object GetFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetInfoString( arg1: FingeprintGenerator32) -> str: 
        """
        GetInfoString( arg1: FingeprintGenerator32) -> str
            Returns a string containing information about the fingerprint generator
            
              RETURNS: an information string
            
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInfoString(RDKit::FingerprintGenerator<unsigned int> const*)
        """
    @staticmethod
    def GetOptions( arg1: FingeprintGenerator32) -> FingerprintOptions: 
        """
        GetOptions( arg1: FingeprintGenerator32) -> FingerprintOptions
            return the fingerprint options object

            C++ signature :
                RDKit::FingerprintArguments* GetOptions(RDKit::FingerprintGenerator<unsigned int>*)
        """
    @staticmethod
    def GetSparseCountFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect: 
        """
        GetSparseCountFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect
            Generates a sparse count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseIntVect containing fingerprint
            
            

            C++ signature :
                RDKit::SparseIntVect<unsigned int>* GetSparseCountFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetSparseFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> SparseBitVect: 
        """
        GetSparseFingerprint( arg1: FingeprintGenerator32, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> SparseBitVect
            Generates a sparse fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseBitVect containing fingerprint
            
            

            C++ signature :
                SparseBitVect* GetSparseFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    pass
class FingeprintGenerator64(Boost.Python.instance):
    @staticmethod
    def GetCountFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect: 
        """
        GetCountFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> UIntSparseIntVect
            Generates a count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseIntVect containing fingerprint
            
            

            C++ signature :
                RDKit::SparseIntVect<unsigned int>* GetCountFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetCountFingerprintAsNumPy( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object: 
        """
        GetCountFingerprintAsNumPy( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object
            Generates a count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a numpy array containing the fingerprint
            
            

            C++ signature :
                boost::python::api::object GetCountFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ExplicitBitVect: 
        """
        GetFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ExplicitBitVect
            Generates a fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a ExplicitBitVect containing fingerprint
            
            

            C++ signature :
                ExplicitBitVect* GetFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetFingerprintAsNumPy( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object: 
        """
        GetFingerprintAsNumPy( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> object
            Generates a fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a numpy array containing the fingerprint
            
            

            C++ signature :
                boost::python::api::object GetFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetInfoString( arg1: FingeprintGenerator64) -> str: 
        """
        GetInfoString( arg1: FingeprintGenerator64) -> str
            Returns a string containing information about the fingerprint generator
            
              RETURNS: an information string
            
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInfoString(RDKit::FingerprintGenerator<unsigned long> const*)
        """
    @staticmethod
    def GetOptions( arg1: FingeprintGenerator64) -> FingerprintOptions: 
        """
        GetOptions( arg1: FingeprintGenerator64) -> FingerprintOptions
            return the fingerprint options object

            C++ signature :
                RDKit::FingerprintArguments* GetOptions(RDKit::FingerprintGenerator<unsigned long>*)
        """
    @staticmethod
    def GetSparseCountFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ULongSparseIntVect: 
        """
        GetSparseCountFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> ULongSparseIntVect
            Generates a sparse count fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseIntVect containing fingerprint
            
            

            C++ signature :
                RDKit::SparseIntVect<unsigned long>* GetSparseCountFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    @staticmethod
    def GetSparseFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> SparseBitVect: 
        """
        GetSparseFingerprint( arg1: FingeprintGenerator64, mol: Mol, fromAtoms: AtomPairsParameters = [], ignoreAtoms: AtomPairsParameters = [], confId: int = -1, customAtomInvariants: AtomPairsParameters = [], customBondInvariants: AtomPairsParameters = [], additionalOutput: AtomPairsParameters = None) -> SparseBitVect
            Generates a sparse fingerprint
            
              ARGUMENTS:
                - mol: molecule to be fingerprinted
                - fromAtoms: indices of atoms to use while generating the fingerprint
                - ignoreAtoms: indices of atoms to exclude while generating the fingerprint
                - confId: 3D confirmation to use, only used by AtomPair fingerprint
                - customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
                - customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
            
                - additionalOutput: AdditionalOutput instance used to return extra information about the bits
            
              RETURNS: a SparseBitVect containing fingerprint
            
            

            C++ signature :
                SparseBitVect* GetSparseFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
    pass
class AtomPairFingerprintOptions(FingerprintOptions, Boost.Python.instance):
    @property
    def maxDistance(self) -> None:
        """
        maximum distance to be included

        :type: None
        """
    @property
    def minDistance(self) -> None:
        """
        minimum distance to be included

        :type: None
        """
    @property
    def use2D(self) -> None:
        """
        use 2D distances

        :type: None
        """
    pass
class MorganFingerprintOptions(FingerprintOptions, Boost.Python.instance):
    @property
    def includeRedundantEnvironments(self) -> None:
        """
        include redundant environments in the fingerprint

        :type: None
        """
    @property
    def onlyNonzeroInvariants(self) -> None:
        """
        use include atoms which have nonzero invariants

        :type: None
        """
    @property
    def radius(self) -> None:
        """
        the radius of the fingerprints to generate

        :type: None
        """
    pass
class RDKitFingerprintOptions(FingerprintOptions, Boost.Python.instance):
    @property
    def branchedPaths(self) -> None:
        """
        generate branched subgraphs, not just linear ones

        :type: None
        """
    @property
    def maxPath(self) -> None:
        """
        maximum path length (in bonds) to be included

        :type: None
        """
    @property
    def minPath(self) -> None:
        """
        minimum path length (in bonds) to be included

        :type: None
        """
    @property
    def useBondOrder(self) -> None:
        """
        include bond orders in the path hashes

        :type: None
        """
    @property
    def useHs(self) -> None:
        """
        use explicit Hs in the paths (if molecule has explicit Hs)

        :type: None
        """
    pass
class TopologicalTorsionFingerprintOptions(FingerprintOptions, Boost.Python.instance):
    @property
    def onlyShortestPaths(self) -> None:
        """
        whether or not to only include paths which are the shortest path between the start and end atoms

        :type: None
        """
    @property
    def torsionAtomCount(self) -> None:
        """
        number of atoms to be included in the paths

        :type: None
        """
    pass
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
def GetCountFPs(molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    GetCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
def GetFPs(molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    GetFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
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
def GetSparseCountFPs(molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    GetSparseCountFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetSparseCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
def GetSparseFPs(molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    GetSparseFPs(molecules: list = [], fpType: FPType = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP) -> list

        C++ signature :
            boost::python::list GetSparseFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
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
AtomPairFP = rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP
MorganFP = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP
RDKitFP = rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP
TopologicalTorsionFP = rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP
