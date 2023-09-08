"""Module containing the core chemistry functionality of the RDKit"""
from __future__ import annotations
import rdkit.Chem.rdchem
import typing
import Boost.Python

__all__ = [
    "ALLOW_CHARGE_SEPARATION",
    "ALLOW_INCOMPLETE_OCTETS",
    "AddMolSubstanceGroup",
    "AllProps",
    "Atom",
    "AtomKekulizeException",
    "AtomMonomerInfo",
    "AtomMonomerType",
    "AtomPDBResidueInfo",
    "AtomProps",
    "AtomSanitizeException",
    "AtomValenceException",
    "Bond",
    "BondDir",
    "BondProps",
    "BondStereo",
    "BondType",
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
    "ChiralType",
    "ClearMolSubstanceGroups",
    "CompositeQueryType",
    "ComputedProps",
    "Conformer",
    "CoordsAsDouble",
    "CreateMolDataSubstanceGroup",
    "CreateMolSubstanceGroup",
    "CreateStereoGroup",
    "EditableMol",
    "FixedMolSizeMolBundle",
    "ForwardStereoGroupIds",
    "GetAtomAlias",
    "GetAtomRLabel",
    "GetAtomValue",
    "GetDefaultPickleProperties",
    "GetMolSubstanceGroupWithIdx",
    "GetMolSubstanceGroups",
    "GetPeriodicTable",
    "GetSupplementalSmilesLabel",
    "HybridizationType",
    "KEKULE_ALL",
    "KekulizeException",
    "Mol",
    "MolBundle",
    "MolBundleCanSerialize",
    "MolProps",
    "MolSanitizeException",
    "NoConformers",
    "NoProps",
    "PeriodicTable",
    "PrivateProps",
    "PropertyPickleOptions",
    "QueryAtom",
    "QueryAtomData",
    "QueryBond",
    "RWMol",
    "ResonanceFlags",
    "ResonanceMolSupplier",
    "ResonanceMolSupplierCallback",
    "RingInfo",
    "STEREO_ABSOLUTE",
    "STEREO_AND",
    "STEREO_OR",
    "SetAtomAlias",
    "SetAtomRLabel",
    "SetAtomValue",
    "SetDefaultPickleProperties",
    "SetSupplementalSmilesLabel",
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
    "UNCONSTRAINED_ANIONS",
    "UNCONSTRAINED_CATIONS",
    "tossit"
]


class Atom(Boost.Python.instance):
    """
    The class to store Atoms.
    Note that, though it is possible to create one, having an Atom on its own
    (i.e not associated with a molecule) is not particularly useful.
    """
    @staticmethod
    def ClearProp( arg1: Atom, arg2: str) -> None: 
        """
        ClearProp( arg1: Atom, arg2: str) -> None
            Removes a particular property from an Atom (does nothing if not already set).
            
              ARGUMENTS:
                - key: the name of the property to be removed.
            

            C++ signature :
                void ClearProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def DescribeQuery( arg1: Atom) -> str: 
        """
        DescribeQuery( arg1: Atom) -> str
            returns a text description of the query. Primarily intended for debugging purposes.
            
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > DescribeQuery(RDKit::Atom const*)
        """
    @staticmethod
    def GetAtomMapNum( arg1: Atom) -> int: 
        """
        GetAtomMapNum( arg1: Atom) -> int
            Gets the atoms map number, returns 0 if not set

            C++ signature :
                int GetAtomMapNum(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetAtomicNum( arg1: Atom) -> int: 
        """
        GetAtomicNum( arg1: Atom) -> int
            Returns the atomic number.

            C++ signature :
                int GetAtomicNum(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetBonds( arg1: Atom) -> tuple: 
        """
        GetBonds( arg1: Atom) -> tuple
            Returns a read-only sequence of the atom's bonds
            

            C++ signature :
                boost::python::tuple GetBonds(RDKit::Atom*)
        """
    @staticmethod
    def GetBoolProp( arg1: Atom, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: Atom, arg2: str) -> bool
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a bool).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                bool GetBoolProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def GetChiralTag( arg1: Atom) -> ChiralType: 
        """
        GetChiralTag( arg1: Atom) -> ChiralType

            C++ signature :
                RDKit::Atom::ChiralType GetChiralTag(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetDegree( arg1: Atom) -> int: 
        """
        GetDegree( arg1: Atom) -> int
            Returns the degree of the atom in the molecule.
            
              The degree of an atom is defined to be its number of
              directly-bonded neighbors.
              The degree is independent of bond orders, but is dependent
                on whether or not Hs are explicit in the graph.
            

            C++ signature :
                unsigned int GetDegree(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetDoubleProp( arg1: Atom, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: Atom, arg2: str) -> float
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                double GetDoubleProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def GetExplicitBitVectProp( arg1: Atom, arg2: str) -> ExplicitBitVect: 
        """
        GetExplicitBitVectProp( arg1: Atom, arg2: str) -> ExplicitBitVect
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a ExplicitBitVect).
            
              RETURNS: an ExplicitBitVect 
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                ExplicitBitVect GetExplicitBitVectProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def GetExplicitValence( arg1: Atom) -> int: 
        """
        GetExplicitValence( arg1: Atom) -> int
            Returns the explicit valence of the atom.
            

            C++ signature :
                int GetExplicitValence(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetFormalCharge( arg1: Atom) -> int: 
        """
        GetFormalCharge( arg1: Atom) -> int

            C++ signature :
                int GetFormalCharge(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetHybridization( arg1: Atom) -> HybridizationType: 
        """
        GetHybridization( arg1: Atom) -> HybridizationType
            Returns the atom's hybridization.
            

            C++ signature :
                RDKit::Atom::HybridizationType GetHybridization(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetIdx( arg1: Atom) -> int: 
        """
        GetIdx( arg1: Atom) -> int
            Returns the atom's index (ordering in the molecule)
            

            C++ signature :
                unsigned int GetIdx(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetImplicitValence( arg1: Atom) -> int: 
        """
        GetImplicitValence( arg1: Atom) -> int
            Returns the number of implicit Hs on the atom.
            

            C++ signature :
                int GetImplicitValence(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetIntProp( arg1: Atom, arg2: str) -> int: 
        """
        GetIntProp( arg1: Atom, arg2: str) -> int
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
              RETURNS: an int
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                int GetIntProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def GetIsAromatic( arg1: Atom) -> bool: 
        """
        GetIsAromatic( arg1: Atom) -> bool

            C++ signature :
                bool GetIsAromatic(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetIsotope( arg1: Atom) -> int: 
        """
        GetIsotope( arg1: Atom) -> int

            C++ signature :
                unsigned int GetIsotope(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetMass( arg1: Atom) -> float: 
        """
        GetMass( arg1: Atom) -> float

            C++ signature :
                double GetMass(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetMonomerInfo( arg1: Atom) -> AtomMonomerInfo: 
        """
        GetMonomerInfo( arg1: Atom) -> AtomMonomerInfo
            Returns the atom's MonomerInfo object, if there is one.
            
            

            C++ signature :
                RDKit::AtomMonomerInfo* GetMonomerInfo(RDKit::Atom*)
        """
    @staticmethod
    def GetNeighbors( arg1: Atom) -> tuple: 
        """
        GetNeighbors( arg1: Atom) -> tuple
            Returns a read-only sequence of the atom's neighbors
            

            C++ signature :
                boost::python::tuple GetNeighbors(RDKit::Atom*)
        """
    @staticmethod
    def GetNoImplicit( arg1: Atom) -> bool: 
        """
        GetNoImplicit( arg1: Atom) -> bool
            Returns whether or not the atom is *allowed* to have implicit Hs.
            

            C++ signature :
                bool GetNoImplicit(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetNumExplicitHs( arg1: Atom) -> int: 
        """
        GetNumExplicitHs( arg1: Atom) -> int

            C++ signature :
                unsigned int GetNumExplicitHs(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetNumImplicitHs( arg1: Atom) -> int: 
        """
        GetNumImplicitHs( arg1: Atom) -> int
            Returns the total number of implicit Hs on the atom.
            

            C++ signature :
                unsigned int GetNumImplicitHs(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetNumRadicalElectrons( arg1: Atom) -> int: 
        """
        GetNumRadicalElectrons( arg1: Atom) -> int

            C++ signature :
                unsigned int GetNumRadicalElectrons(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetOwningMol( arg1: Atom) -> Mol: 
        """
        GetOwningMol( arg1: Atom) -> Mol
            Returns the Mol that owns this atom.
            

            C++ signature :
                RDKit::ROMol {lvalue} GetOwningMol(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetPDBResidueInfo( arg1: Atom) -> AtomPDBResidueInfo: 
        """
        GetPDBResidueInfo( arg1: Atom) -> AtomPDBResidueInfo
            Returns the atom's MonomerInfo object, if there is one.
            
            

            C++ signature :
                RDKit::AtomPDBResidueInfo* GetPDBResidueInfo(RDKit::Atom*)
        """
    def GetProp(self, key: str, autoConvert: bool = False) -> object: 
        """
        GetProp( self: Atom, key: str, autoConvert: bool = False) -> object
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                boost::python::api::object GetProp(RDKit::Atom const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: Atom, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a list of the properties set on the Atom.
            
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Atom {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: Atom, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict
            Returns a dictionary of the properties set on the Atom.
             n.b. some properties cannot be converted to python types.
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::Atom [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    def GetQueryType( arg1: Atom) -> str: 
        """
        GetQueryType( arg1: Atom) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetQueryType(RDKit::Atom {lvalue})
        """
    def GetSmarts(self, doKekule: bool = False, allHsExplicit: bool = False, isomericSmiles: bool = True) -> str: 
        """
        GetSmarts( self: Atom, doKekule: bool = False, allHsExplicit: bool = False, isomericSmiles: bool = True) -> str
            returns the SMARTS (or SMILES) string for an Atom
            
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSmarts(RDKit::Atom const* [,bool=False [,bool=False [,bool=True]]])
        """
    @staticmethod
    def GetSymbol( arg1: Atom) -> str: 
        """
        GetSymbol( arg1: Atom) -> str
            Returns the atomic symbol (a string)
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSymbol(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetTotalDegree( arg1: Atom) -> int: 
        """
        GetTotalDegree( arg1: Atom) -> int
            Returns the degree of the atom in the molecule including Hs.
            
              The degree of an atom is defined to be its number of
              directly-bonded neighbors.
              The degree is independent of bond orders.
            

            C++ signature :
                unsigned int GetTotalDegree(RDKit::Atom {lvalue})
        """
    def GetTotalNumHs(self, includeNeighbors: bool = False) -> int: 
        """
        GetTotalNumHs( self: Atom, includeNeighbors: bool = False) -> int
            Returns the total number of Hs (explicit and implicit) on the atom.
            
              ARGUMENTS:
            
                - includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.
                  Defaults to 0.
            

            C++ signature :
                unsigned int GetTotalNumHs(RDKit::Atom {lvalue} [,bool=False])
        """
    @staticmethod
    def GetTotalValence( arg1: Atom) -> int: 
        """
        GetTotalValence( arg1: Atom) -> int
            Returns the total valence (explicit + implicit) of the atom.
            
            

            C++ signature :
                unsigned int GetTotalValence(RDKit::Atom {lvalue})
        """
    @staticmethod
    def GetUnsignedProp( arg1: Atom, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: Atom, arg2: str) -> int
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
              RETURNS: an integer (Python has no unsigned type)
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def HasOwningMol( arg1: Atom) -> bool: 
        """
        HasOwningMol( arg1: Atom) -> bool
            Returns whether or not this instance belongs to a molecule.
            

            C++ signature :
                bool HasOwningMol(RDKit::Atom {lvalue})
        """
    @staticmethod
    def HasProp( arg1: Atom, arg2: str) -> int: 
        """
        HasProp( arg1: Atom, arg2: str) -> int
            Queries a Atom to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            

            C++ signature :
                int HasProp(RDKit::Atom const*,char const*)
        """
    @staticmethod
    def HasQuery( arg1: Atom) -> bool: 
        """
        HasQuery( arg1: Atom) -> bool
            Returns whether or not the atom has an associated query
            
            

            C++ signature :
                bool HasQuery(RDKit::Atom {lvalue})
        """
    @staticmethod
    def InvertChirality( arg1: Atom) -> bool: 
        """
        InvertChirality( arg1: Atom) -> bool

            C++ signature :
                bool InvertChirality(RDKit::Atom {lvalue})
        """
    @staticmethod
    def IsInRing( arg1: Atom) -> bool: 
        """
        IsInRing( arg1: Atom) -> bool
            Returns whether or not the atom is in a ring
            
            

            C++ signature :
                bool IsInRing(RDKit::Atom const*)
        """
    @staticmethod
    def IsInRingSize( arg1: Atom, arg2: int) -> bool: 
        """
        IsInRingSize( arg1: Atom, arg2: int) -> bool
            Returns whether or not the atom is in a ring of a particular size.
            
              ARGUMENTS:
                - size: the ring size to look for
            

            C++ signature :
                bool IsInRingSize(RDKit::Atom const*,int)
        """
    @staticmethod
    def Match( arg1: Atom, arg2: Atom) -> bool: 
        """
        Match( arg1: Atom, arg2: Atom) -> bool
            Returns whether or not this atom matches another Atom.
            
              Each Atom (or query Atom) has a query function which is
              used for this type of matching.
            
              ARGUMENTS:
                - other: the other Atom to which to compare
            

            C++ signature :
                bool Match(RDKit::Atom {lvalue},RDKit::Atom const*)
        """
    def NeedsUpdatePropertyCache(self) -> bool: 
        """
        NeedsUpdatePropertyCache( self: Atom) -> bool
            Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
            
            

            C++ signature :
                bool NeedsUpdatePropertyCache(RDKit::Atom {lvalue})
        """
    def SetAtomMapNum(self, mapno: int, strict: bool = False) -> None: 
        """
        SetAtomMapNum( self: Atom, mapno: int, strict: bool = False) -> None
            Sets the atoms map number, a value of 0 clears the atom map

            C++ signature :
                void SetAtomMapNum(RDKit::Atom {lvalue},int [,bool=False])
        """
    @staticmethod
    def SetAtomicNum( arg1: Atom, arg2: int) -> None: 
        """
        SetAtomicNum( arg1: Atom, arg2: int) -> None
            Sets the atomic number, takes an integer value as an argument

            C++ signature :
                void SetAtomicNum(RDKit::Atom {lvalue},int)
        """
    def SetBoolProp(self, key: str, val: bool) -> None: 
        """
        SetBoolProp( self: Atom, key: str, val: bool) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a bool).
                - value: the property value (a bool).
            
            

            C++ signature :
                void SetBoolProp(RDKit::Atom const*,char const*,bool)
        """
    @staticmethod
    def SetChiralTag( arg1: Atom, arg2: ChiralType) -> None: 
        """
        SetChiralTag( arg1: Atom, arg2: ChiralType) -> None

            C++ signature :
                void SetChiralTag(RDKit::Atom {lvalue},RDKit::Atom::ChiralType)
        """
    def SetDoubleProp(self, key: str, val: float) -> None: 
        """
        SetDoubleProp( self: Atom, key: str, val: float) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a double).
                - value: the property value (a double).
            
            

            C++ signature :
                void SetDoubleProp(RDKit::Atom const*,char const*,double)
        """
    def SetExplicitBitVectProp(self, key: str, val: ExplicitBitVect) -> None: 
        """
        SetExplicitBitVectProp( self: Atom, key: str, val: ExplicitBitVect) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (an ExplicitBitVect).
                - value: the property value (an ExplicitBitVect).
            
            

            C++ signature :
                void SetExplicitBitVectProp(RDKit::Atom const*,char const*,ExplicitBitVect)
        """
    @staticmethod
    def SetFormalCharge( arg1: Atom, arg2: int) -> None: 
        """
        SetFormalCharge( arg1: Atom, arg2: int) -> None

            C++ signature :
                void SetFormalCharge(RDKit::Atom {lvalue},int)
        """
    @staticmethod
    def SetHybridization( arg1: Atom, arg2: HybridizationType) -> None: 
        """
        SetHybridization( arg1: Atom, arg2: HybridizationType) -> None
            Sets the hybridization of the atom.
              The argument should be a HybridizationType
            

            C++ signature :
                void SetHybridization(RDKit::Atom {lvalue},RDKit::Atom::HybridizationType)
        """
    def SetIntProp(self, key: str, val: int) -> None: 
        """
        SetIntProp( self: Atom, key: str, val: int) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a int).
                - value: the property value (a int).
            
            

            C++ signature :
                void SetIntProp(RDKit::Atom const*,char const*,int)
        """
    @staticmethod
    def SetIsAromatic( arg1: Atom, arg2: bool) -> None: 
        """
        SetIsAromatic( arg1: Atom, arg2: bool) -> None

            C++ signature :
                void SetIsAromatic(RDKit::Atom {lvalue},bool)
        """
    @staticmethod
    def SetIsotope( arg1: Atom, arg2: int) -> None: 
        """
        SetIsotope( arg1: Atom, arg2: int) -> None

            C++ signature :
                void SetIsotope(RDKit::Atom {lvalue},unsigned int)
        """
    @staticmethod
    def SetMonomerInfo( arg1: Atom, arg2: AtomMonomerInfo) -> None: 
        """
        SetMonomerInfo( arg1: Atom, arg2: AtomMonomerInfo) -> None
            Sets the atom's MonomerInfo object.
            
            

            C++ signature :
                void SetMonomerInfo(RDKit::Atom*,RDKit::AtomMonomerInfo const*)
        """
    @staticmethod
    def SetNoImplicit( arg1: Atom, arg2: bool) -> None: 
        """
        SetNoImplicit( arg1: Atom, arg2: bool) -> None
            Sets a marker on the atom that *disallows* implicit Hs.
              This holds even if the atom would otherwise have implicit Hs added.
            

            C++ signature :
                void SetNoImplicit(RDKit::Atom {lvalue},bool)
        """
    @staticmethod
    def SetNumExplicitHs( arg1: Atom, arg2: int) -> None: 
        """
        SetNumExplicitHs( arg1: Atom, arg2: int) -> None

            C++ signature :
                void SetNumExplicitHs(RDKit::Atom {lvalue},unsigned int)
        """
    @staticmethod
    def SetNumRadicalElectrons( arg1: Atom, arg2: int) -> None: 
        """
        SetNumRadicalElectrons( arg1: Atom, arg2: int) -> None

            C++ signature :
                void SetNumRadicalElectrons(RDKit::Atom {lvalue},unsigned int)
        """
    @staticmethod
    def SetPDBResidueInfo( arg1: Atom, arg2: AtomMonomerInfo) -> None: 
        """
        SetPDBResidueInfo( arg1: Atom, arg2: AtomMonomerInfo) -> None
            Sets the atom's MonomerInfo object.
            
            

            C++ signature :
                void SetPDBResidueInfo(RDKit::Atom*,RDKit::AtomMonomerInfo const*)
        """
    def SetProp(self, key: str, val: str) -> None: 
        """
        SetProp( self: Atom, key: str, val: str) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
            
            

            C++ signature :
                void SetProp(RDKit::Atom const*,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    def SetUnsignedProp(self, key: str, val: int) -> None: 
        """
        SetUnsignedProp( self: Atom, key: str, val: int) -> None
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned integer).
                - value: the property value (a int >= 0).
            
            

            C++ signature :
                void SetUnsignedProp(RDKit::Atom const*,char const*,unsigned int)
        """
    def UpdatePropertyCache(self, strict: bool = True) -> None: 
        """
        UpdatePropertyCache( self: Atom, strict: bool = True) -> None
            Regenerates computed properties like implicit valence and ring information.
            
            

            C++ signature :
                void UpdatePropertyCache(RDKit::Atom {lvalue} [,bool=True])
        """
    @staticmethod
    def __copy__( arg1: Atom) -> Atom: 
        """
        __copy__( arg1: Atom) -> Atom
            Create a copy of the atom

            C++ signature :
                RDKit::Atom* __copy__(RDKit::Atom {lvalue})
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: 
        """
        __init__( arg1: object, arg2: str) -> None

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,RDKit::Atom)

            C++ signature :
                void __init__(_object*,unsigned int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: Atom) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: ...
    __instance_size__ = 96
    pass
class MolSanitizeException(ValueError, Exception, BaseException):
    pass
class AtomMonomerInfo(Boost.Python.instance):
    """
    The class to store monomer information attached to Atoms
    """
    @staticmethod
    def GetMonomerType( arg1: AtomMonomerInfo) -> AtomMonomerType: 
        """
        GetMonomerType( arg1: AtomMonomerInfo) -> AtomMonomerType

            C++ signature :
                RDKit::AtomMonomerInfo::AtomMonomerType GetMonomerType(RDKit::AtomMonomerInfo {lvalue})
        """
    @staticmethod
    def GetName( arg1: AtomMonomerInfo) -> str: 
        """
        GetName( arg1: AtomMonomerInfo) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetName(RDKit::AtomMonomerInfo {lvalue})
        """
    @staticmethod
    def SetMonomerType( arg1: AtomMonomerInfo, arg2: AtomMonomerType) -> None: 
        """
        SetMonomerType( arg1: AtomMonomerInfo, arg2: AtomMonomerType) -> None

            C++ signature :
                void SetMonomerType(RDKit::AtomMonomerInfo {lvalue},RDKit::AtomMonomerInfo::AtomMonomerType)
        """
    @staticmethod
    def SetName( arg1: AtomMonomerInfo, arg2: str) -> None: 
        """
        SetName( arg1: AtomMonomerInfo, arg2: str) -> None

            C++ signature :
                void SetName(RDKit::AtomMonomerInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,RDKit::AtomMonomerInfo::AtomMonomerType [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, type: AtomMonomerType, name: str = '') -> None: ...
    __instance_size__ = 72
    pass
class AtomMonomerType(Boost.Python.enum, int):
    OTHER = rdkit.Chem.rdchem.AtomMonomerType.OTHER
    PDBRESIDUE = rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE
    UNKNOWN = rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN
    __slots__ = ()
    names = {'UNKNOWN': rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN, 'PDBRESIDUE': rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE, 'OTHER': rdkit.Chem.rdchem.AtomMonomerType.OTHER}
    values = {0: rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN, 1: rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE, 2: rdkit.Chem.rdchem.AtomMonomerType.OTHER}
    pass
class AtomPDBResidueInfo(AtomMonomerInfo, Boost.Python.instance):
    """
    The class to store PDB residue information attached to Atoms
    """
    @staticmethod
    def GetAltLoc( arg1: AtomPDBResidueInfo) -> str: 
        """
        GetAltLoc( arg1: AtomPDBResidueInfo) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetAltLoc(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetChainId( arg1: AtomPDBResidueInfo) -> str: 
        """
        GetChainId( arg1: AtomPDBResidueInfo) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetChainId(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetInsertionCode( arg1: AtomPDBResidueInfo) -> str: 
        """
        GetInsertionCode( arg1: AtomPDBResidueInfo) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInsertionCode(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetIsHeteroAtom( arg1: AtomPDBResidueInfo) -> bool: 
        """
        GetIsHeteroAtom( arg1: AtomPDBResidueInfo) -> bool

            C++ signature :
                bool GetIsHeteroAtom(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetOccupancy( arg1: AtomPDBResidueInfo) -> float: 
        """
        GetOccupancy( arg1: AtomPDBResidueInfo) -> float

            C++ signature :
                double GetOccupancy(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetResidueName( arg1: AtomPDBResidueInfo) -> str: 
        """
        GetResidueName( arg1: AtomPDBResidueInfo) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetResidueName(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetResidueNumber( arg1: AtomPDBResidueInfo) -> int: 
        """
        GetResidueNumber( arg1: AtomPDBResidueInfo) -> int

            C++ signature :
                int GetResidueNumber(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetSecondaryStructure( arg1: AtomPDBResidueInfo) -> int: 
        """
        GetSecondaryStructure( arg1: AtomPDBResidueInfo) -> int

            C++ signature :
                unsigned int GetSecondaryStructure(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetSegmentNumber( arg1: AtomPDBResidueInfo) -> int: 
        """
        GetSegmentNumber( arg1: AtomPDBResidueInfo) -> int

            C++ signature :
                unsigned int GetSegmentNumber(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetSerialNumber( arg1: AtomPDBResidueInfo) -> int: 
        """
        GetSerialNumber( arg1: AtomPDBResidueInfo) -> int

            C++ signature :
                int GetSerialNumber(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def GetTempFactor( arg1: AtomPDBResidueInfo) -> float: 
        """
        GetTempFactor( arg1: AtomPDBResidueInfo) -> float

            C++ signature :
                double GetTempFactor(RDKit::AtomPDBResidueInfo {lvalue})
        """
    @staticmethod
    def SetAltLoc( arg1: AtomPDBResidueInfo, arg2: str) -> None: 
        """
        SetAltLoc( arg1: AtomPDBResidueInfo, arg2: str) -> None

            C++ signature :
                void SetAltLoc(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetChainId( arg1: AtomPDBResidueInfo, arg2: str) -> None: 
        """
        SetChainId( arg1: AtomPDBResidueInfo, arg2: str) -> None

            C++ signature :
                void SetChainId(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetInsertionCode( arg1: AtomPDBResidueInfo, arg2: str) -> None: 
        """
        SetInsertionCode( arg1: AtomPDBResidueInfo, arg2: str) -> None

            C++ signature :
                void SetInsertionCode(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetIsHeteroAtom( arg1: AtomPDBResidueInfo, arg2: bool) -> None: 
        """
        SetIsHeteroAtom( arg1: AtomPDBResidueInfo, arg2: bool) -> None

            C++ signature :
                void SetIsHeteroAtom(RDKit::AtomPDBResidueInfo {lvalue},bool)
        """
    @staticmethod
    def SetOccupancy( arg1: AtomPDBResidueInfo, arg2: float) -> None: 
        """
        SetOccupancy( arg1: AtomPDBResidueInfo, arg2: float) -> None

            C++ signature :
                void SetOccupancy(RDKit::AtomPDBResidueInfo {lvalue},double)
        """
    @staticmethod
    def SetResidueName( arg1: AtomPDBResidueInfo, arg2: str) -> None: 
        """
        SetResidueName( arg1: AtomPDBResidueInfo, arg2: str) -> None

            C++ signature :
                void SetResidueName(RDKit::AtomPDBResidueInfo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetResidueNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None: 
        """
        SetResidueNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None

            C++ signature :
                void SetResidueNumber(RDKit::AtomPDBResidueInfo {lvalue},int)
        """
    @staticmethod
    def SetSecondaryStructure( arg1: AtomPDBResidueInfo, arg2: int) -> None: 
        """
        SetSecondaryStructure( arg1: AtomPDBResidueInfo, arg2: int) -> None

            C++ signature :
                void SetSecondaryStructure(RDKit::AtomPDBResidueInfo {lvalue},unsigned int)
        """
    @staticmethod
    def SetSegmentNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None: 
        """
        SetSegmentNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None

            C++ signature :
                void SetSegmentNumber(RDKit::AtomPDBResidueInfo {lvalue},unsigned int)
        """
    @staticmethod
    def SetSerialNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None: 
        """
        SetSerialNumber( arg1: AtomPDBResidueInfo, arg2: int) -> None

            C++ signature :
                void SetSerialNumber(RDKit::AtomPDBResidueInfo {lvalue},int)
        """
    @staticmethod
    def SetTempFactor( arg1: AtomPDBResidueInfo, arg2: float) -> None: 
        """
        SetTempFactor( arg1: AtomPDBResidueInfo, arg2: float) -> None

            C++ signature :
                void SetTempFactor(RDKit::AtomPDBResidueInfo {lvalue},double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,int=1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,int=0 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,double=1.0 [,double=0.0 [,bool=False [,unsigned int=0 [,unsigned int=0]]]]]]]]]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, atomName: str, serialNumber: int = 1, altLoc: str = '', residueName: str = '', residueNumber: int = 0, chainId: str = '', insertionCode: str = '', occupancy: float = 1.0, tempFactor: float = 0.0, isHeteroAtom: bool = False, secondaryStructure: int = 0, segmentNumber: int = 0) -> None: ...
    __instance_size__ = 248
    pass
class AtomSanitizeException(MolSanitizeException, ValueError, Exception, BaseException):
    pass
class AtomValenceException(AtomSanitizeException, MolSanitizeException, ValueError, Exception, BaseException):
    pass
class Bond(Boost.Python.instance):
    """
    The class to store Bonds.
    Note: unlike Atoms, is it currently impossible to construct Bonds from
    Python.
    """
    @staticmethod
    def ClearProp( arg1: Bond, arg2: str) -> None: 
        """
        ClearProp( arg1: Bond, arg2: str) -> None
            Removes a particular property from an Bond (does nothing if not already set).
            
              ARGUMENTS:
                - key: the name of the property to be removed.
            

            C++ signature :
                void ClearProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def DescribeQuery( arg1: Bond) -> str: 
        """
        DescribeQuery( arg1: Bond) -> str
            returns a text description of the query. Primarily intended for debugging purposes.
            
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > DescribeQuery(RDKit::Bond const*)
        """
    @staticmethod
    def GetBeginAtom( arg1: Bond) -> Atom: 
        """
        GetBeginAtom( arg1: Bond) -> Atom
            Returns the bond's first atom.
            

            C++ signature :
                RDKit::Atom* GetBeginAtom(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetBeginAtomIdx( arg1: Bond) -> int: 
        """
        GetBeginAtomIdx( arg1: Bond) -> int
            Returns the index of the bond's first atom.
            

            C++ signature :
                unsigned int GetBeginAtomIdx(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetBondDir( arg1: Bond) -> BondDir: 
        """
        GetBondDir( arg1: Bond) -> BondDir
            Returns the type of the bond as a BondDir
            

            C++ signature :
                RDKit::Bond::BondDir GetBondDir(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetBondType( arg1: Bond) -> BondType: 
        """
        GetBondType( arg1: Bond) -> BondType
            Returns the type of the bond as a BondType
            

            C++ signature :
                RDKit::Bond::BondType GetBondType(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetBondTypeAsDouble( arg1: Bond) -> float: 
        """
        GetBondTypeAsDouble( arg1: Bond) -> float
            Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
            

            C++ signature :
                double GetBondTypeAsDouble(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetBoolProp( arg1: Bond, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: Bond, arg2: str) -> bool
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a boolean).
            
              RETURNS: a boolean
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                bool GetBoolProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def GetDoubleProp( arg1: Bond, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: Bond, arg2: str) -> float
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                double GetDoubleProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def GetEndAtom( arg1: Bond) -> Atom: 
        """
        GetEndAtom( arg1: Bond) -> Atom
            Returns the bond's second atom.
            

            C++ signature :
                RDKit::Atom* GetEndAtom(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetEndAtomIdx( arg1: Bond) -> int: 
        """
        GetEndAtomIdx( arg1: Bond) -> int
            Returns the index of the bond's first atom.
            

            C++ signature :
                unsigned int GetEndAtomIdx(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetIdx( arg1: Bond) -> int: 
        """
        GetIdx( arg1: Bond) -> int
            Returns the bond's index (ordering in the molecule)
            

            C++ signature :
                unsigned int GetIdx(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetIntProp( arg1: Bond, arg2: str) -> int: 
        """
        GetIntProp( arg1: Bond, arg2: str) -> int
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
              RETURNS: an int
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                int GetIntProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def GetIsAromatic( arg1: Bond) -> bool: 
        """
        GetIsAromatic( arg1: Bond) -> bool

            C++ signature :
                bool GetIsAromatic(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetIsConjugated( arg1: Bond) -> bool: 
        """
        GetIsConjugated( arg1: Bond) -> bool
            Returns whether or not the bond is considered to be conjugated.

            C++ signature :
                bool GetIsConjugated(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetOtherAtom( arg1: Bond, arg2: Atom) -> Atom: 
        """
        GetOtherAtom( arg1: Bond, arg2: Atom) -> Atom
            Given one of the bond's atoms, returns the other one.
            

            C++ signature :
                RDKit::Atom* GetOtherAtom(RDKit::Bond {lvalue},RDKit::Atom const*)
        """
    @staticmethod
    def GetOtherAtomIdx( arg1: Bond, arg2: int) -> int: 
        """
        GetOtherAtomIdx( arg1: Bond, arg2: int) -> int
            Given the index of one of the bond's atoms, returns the
            index of the other.
            

            C++ signature :
                unsigned int GetOtherAtomIdx(RDKit::Bond {lvalue},unsigned int)
        """
    @staticmethod
    def GetOwningMol( arg1: Bond) -> Mol: 
        """
        GetOwningMol( arg1: Bond) -> Mol
            Returns the Mol that owns this bond.
            

            C++ signature :
                RDKit::ROMol {lvalue} GetOwningMol(RDKit::Bond {lvalue})
        """
    def GetProp(self, key: str, autoConvert: bool = False) -> object: 
        """
        GetProp( self: Bond, key: str, autoConvert: bool = False) -> object
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                boost::python::api::object GetProp(RDKit::Bond const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: Bond, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a list of the properties set on the Bond.
            
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Bond {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: Bond, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict
            Returns a dictionary of the properties set on the Bond.
             n.b. some properties cannot be converted to python types.
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::Bond [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    def GetSmarts( bond: Bond, allBondsExplicit: bool = False) -> str: 
        """
        GetSmarts( bond: Bond, allBondsExplicit: bool = False) -> str
            returns the SMARTS (or SMILES) string for a Bond

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSmarts(RDKit::Bond const* [,bool=False])
        """
    @staticmethod
    def GetStereo( arg1: Bond) -> BondStereo: 
        """
        GetStereo( arg1: Bond) -> BondStereo
            Returns the stereo configuration of the bond as a BondStereo
            

            C++ signature :
                RDKit::Bond::BondStereo GetStereo(RDKit::Bond {lvalue})
        """
    @staticmethod
    def GetStereoAtoms( arg1: Bond) -> _vecti: 
        """
        GetStereoAtoms( arg1: Bond) -> _vecti
            Returns the indices of the atoms setting this bond's stereochemistry.
            

            C++ signature :
                std::vector<int, std::allocator<int> > GetStereoAtoms(RDKit::Bond const*)
        """
    @staticmethod
    def GetUnsignedProp( arg1: Bond, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: Bond, arg2: str) -> int
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
              RETURNS: an int (Python has no unsigned type)
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def GetValenceContrib( arg1: Bond, arg2: Atom) -> float: 
        """
        GetValenceContrib( arg1: Bond, arg2: Atom) -> float
            Returns the contribution of the bond to the valence of an Atom.
            
              ARGUMENTS:
            
                - atom: the Atom to consider.
            

            C++ signature :
                double GetValenceContrib(RDKit::Bond {lvalue},RDKit::Atom const*)
        """
    @staticmethod
    def HasOwningMol( arg1: Bond) -> bool: 
        """
        HasOwningMol( arg1: Bond) -> bool
            Returns whether or not this instance belongs to a molecule.
            

            C++ signature :
                bool HasOwningMol(RDKit::Bond {lvalue})
        """
    @staticmethod
    def HasProp( arg1: Bond, arg2: str) -> int: 
        """
        HasProp( arg1: Bond, arg2: str) -> int
            Queries a Bond to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            

            C++ signature :
                int HasProp(RDKit::Bond const*,char const*)
        """
    @staticmethod
    def HasQuery( arg1: Bond) -> bool: 
        """
        HasQuery( arg1: Bond) -> bool
            Returns whether or not the bond has an associated query
            
            

            C++ signature :
                bool HasQuery(RDKit::Bond {lvalue})
        """
    @staticmethod
    def IsInRing( arg1: Bond) -> bool: 
        """
        IsInRing( arg1: Bond) -> bool
            Returns whether or not the bond is in a ring of any size.
            
            

            C++ signature :
                bool IsInRing(RDKit::Bond const*)
        """
    @staticmethod
    def IsInRingSize( arg1: Bond, arg2: int) -> bool: 
        """
        IsInRingSize( arg1: Bond, arg2: int) -> bool
            Returns whether or not the bond is in a ring of a particular size.
            
              ARGUMENTS:
                - size: the ring size to look for
            

            C++ signature :
                bool IsInRingSize(RDKit::Bond const*,int)
        """
    @staticmethod
    def Match( arg1: Bond, arg2: Bond) -> bool: 
        """
        Match( arg1: Bond, arg2: Bond) -> bool
            Returns whether or not this bond matches another Bond.
            
              Each Bond (or query Bond) has a query function which is
              used for this type of matching.
            
              ARGUMENTS:
                - other: the other Bond to which to compare
            

            C++ signature :
                bool Match(RDKit::Bond {lvalue},RDKit::Bond const*)
        """
    @staticmethod
    def SetBondDir( arg1: Bond, arg2: BondDir) -> None: 
        """
        SetBondDir( arg1: Bond, arg2: BondDir) -> None
            Set the type of the bond as a BondDir
            

            C++ signature :
                void SetBondDir(RDKit::Bond {lvalue},RDKit::Bond::BondDir)
        """
    @staticmethod
    def SetBondType( arg1: Bond, arg2: BondType) -> None: 
        """
        SetBondType( arg1: Bond, arg2: BondType) -> None
            Set the type of the bond as a BondType
            

            C++ signature :
                void SetBondType(RDKit::Bond {lvalue},RDKit::Bond::BondType)
        """
    def SetBoolProp(self, key: str, val: bool) -> None: 
        """
        SetBoolProp( self: Bond, key: str, val: bool) -> None
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a boolean).
            
            

            C++ signature :
                void SetBoolProp(RDKit::Bond const*,char const*,bool)
        """
    def SetDoubleProp(self, key: str, val: float) -> None: 
        """
        SetDoubleProp( self: Bond, key: str, val: float) -> None
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a double).
            
            

            C++ signature :
                void SetDoubleProp(RDKit::Bond const*,char const*,double)
        """
    def SetIntProp(self, key: str, val: int) -> None: 
        """
        SetIntProp( self: Bond, key: str, val: int) -> None
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (an int).
            
            

            C++ signature :
                void SetIntProp(RDKit::Bond const*,char const*,int)
        """
    @staticmethod
    def SetIsAromatic( arg1: Bond, arg2: bool) -> None: 
        """
        SetIsAromatic( arg1: Bond, arg2: bool) -> None

            C++ signature :
                void SetIsAromatic(RDKit::Bond {lvalue},bool)
        """
    @staticmethod
    def SetIsConjugated( arg1: Bond, arg2: bool) -> None: 
        """
        SetIsConjugated( arg1: Bond, arg2: bool) -> None

            C++ signature :
                void SetIsConjugated(RDKit::Bond {lvalue},bool)
        """
    def SetProp(self, key: str, val: str) -> None: 
        """
        SetProp( self: Bond, key: str, val: str) -> None
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
            
            

            C++ signature :
                void SetProp(RDKit::Bond const*,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetStereo( arg1: Bond, arg2: BondStereo) -> None: 
        """
        SetStereo( arg1: Bond, arg2: BondStereo) -> None
            Set the stereo configuration of the bond as a BondStereo
            

            C++ signature :
                void SetStereo(RDKit::Bond {lvalue},RDKit::Bond::BondStereo)
        """
    @staticmethod
    def SetStereoAtoms( arg1: Bond, arg2: int, arg3: int) -> None: 
        """
        SetStereoAtoms( arg1: Bond, arg2: int, arg3: int) -> None
            Set the indices of the atoms setting this bond's stereochemistry.
            

            C++ signature :
                void SetStereoAtoms(RDKit::Bond {lvalue},unsigned int,unsigned int)
        """
    def SetUnsignedProp(self, key: str, val: int) -> None: 
        """
        SetUnsignedProp( self: Bond, key: str, val: int) -> None
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (an int >= 0).
            
            

            C++ signature :
                void SetUnsignedProp(RDKit::Bond const*,char const*,unsigned int)
        """
    pass
class BondDir(Boost.Python.enum, int):
    BEGINDASH = rdkit.Chem.rdchem.BondDir.BEGINDASH
    BEGINWEDGE = rdkit.Chem.rdchem.BondDir.BEGINWEDGE
    EITHERDOUBLE = rdkit.Chem.rdchem.BondDir.EITHERDOUBLE
    ENDDOWNRIGHT = rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT
    ENDUPRIGHT = rdkit.Chem.rdchem.BondDir.ENDUPRIGHT
    NONE = rdkit.Chem.rdchem.BondDir.NONE
    UNKNOWN = rdkit.Chem.rdchem.BondDir.UNKNOWN
    __slots__ = ()
    names = {'NONE': rdkit.Chem.rdchem.BondDir.NONE, 'BEGINWEDGE': rdkit.Chem.rdchem.BondDir.BEGINWEDGE, 'BEGINDASH': rdkit.Chem.rdchem.BondDir.BEGINDASH, 'ENDDOWNRIGHT': rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT, 'ENDUPRIGHT': rdkit.Chem.rdchem.BondDir.ENDUPRIGHT, 'EITHERDOUBLE': rdkit.Chem.rdchem.BondDir.EITHERDOUBLE, 'UNKNOWN': rdkit.Chem.rdchem.BondDir.UNKNOWN}
    values = {0: rdkit.Chem.rdchem.BondDir.NONE, 1: rdkit.Chem.rdchem.BondDir.BEGINWEDGE, 2: rdkit.Chem.rdchem.BondDir.BEGINDASH, 3: rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT, 4: rdkit.Chem.rdchem.BondDir.ENDUPRIGHT, 5: rdkit.Chem.rdchem.BondDir.EITHERDOUBLE, 6: rdkit.Chem.rdchem.BondDir.UNKNOWN}
    pass
class BondStereo(Boost.Python.enum, int):
    STEREOANY = rdkit.Chem.rdchem.BondStereo.STEREOANY
    STEREOCIS = rdkit.Chem.rdchem.BondStereo.STEREOCIS
    STEREOE = rdkit.Chem.rdchem.BondStereo.STEREOE
    STEREONONE = rdkit.Chem.rdchem.BondStereo.STEREONONE
    STEREOTRANS = rdkit.Chem.rdchem.BondStereo.STEREOTRANS
    STEREOZ = rdkit.Chem.rdchem.BondStereo.STEREOZ
    __slots__ = ()
    names = {'STEREONONE': rdkit.Chem.rdchem.BondStereo.STEREONONE, 'STEREOANY': rdkit.Chem.rdchem.BondStereo.STEREOANY, 'STEREOZ': rdkit.Chem.rdchem.BondStereo.STEREOZ, 'STEREOE': rdkit.Chem.rdchem.BondStereo.STEREOE, 'STEREOCIS': rdkit.Chem.rdchem.BondStereo.STEREOCIS, 'STEREOTRANS': rdkit.Chem.rdchem.BondStereo.STEREOTRANS}
    values = {0: rdkit.Chem.rdchem.BondStereo.STEREONONE, 1: rdkit.Chem.rdchem.BondStereo.STEREOANY, 2: rdkit.Chem.rdchem.BondStereo.STEREOZ, 3: rdkit.Chem.rdchem.BondStereo.STEREOE, 4: rdkit.Chem.rdchem.BondStereo.STEREOCIS, 5: rdkit.Chem.rdchem.BondStereo.STEREOTRANS}
    pass
class BondType(Boost.Python.enum, int):
    AROMATIC = rdkit.Chem.rdchem.BondType.AROMATIC
    DATIVE = rdkit.Chem.rdchem.BondType.DATIVE
    DATIVEL = rdkit.Chem.rdchem.BondType.DATIVEL
    DATIVEONE = rdkit.Chem.rdchem.BondType.DATIVEONE
    DATIVER = rdkit.Chem.rdchem.BondType.DATIVER
    DOUBLE = rdkit.Chem.rdchem.BondType.DOUBLE
    FIVEANDAHALF = rdkit.Chem.rdchem.BondType.FIVEANDAHALF
    FOURANDAHALF = rdkit.Chem.rdchem.BondType.FOURANDAHALF
    HEXTUPLE = rdkit.Chem.rdchem.BondType.HEXTUPLE
    HYDROGEN = rdkit.Chem.rdchem.BondType.HYDROGEN
    IONIC = rdkit.Chem.rdchem.BondType.IONIC
    ONEANDAHALF = rdkit.Chem.rdchem.BondType.ONEANDAHALF
    OTHER = rdkit.Chem.rdchem.BondType.OTHER
    QUADRUPLE = rdkit.Chem.rdchem.BondType.QUADRUPLE
    QUINTUPLE = rdkit.Chem.rdchem.BondType.QUINTUPLE
    SINGLE = rdkit.Chem.rdchem.BondType.SINGLE
    THREEANDAHALF = rdkit.Chem.rdchem.BondType.THREEANDAHALF
    THREECENTER = rdkit.Chem.rdchem.BondType.THREECENTER
    TRIPLE = rdkit.Chem.rdchem.BondType.TRIPLE
    TWOANDAHALF = rdkit.Chem.rdchem.BondType.TWOANDAHALF
    UNSPECIFIED = rdkit.Chem.rdchem.BondType.UNSPECIFIED
    ZERO = rdkit.Chem.rdchem.BondType.ZERO
    __slots__ = ()
    names = {'UNSPECIFIED': rdkit.Chem.rdchem.BondType.UNSPECIFIED, 'SINGLE': rdkit.Chem.rdchem.BondType.SINGLE, 'DOUBLE': rdkit.Chem.rdchem.BondType.DOUBLE, 'TRIPLE': rdkit.Chem.rdchem.BondType.TRIPLE, 'QUADRUPLE': rdkit.Chem.rdchem.BondType.QUADRUPLE, 'QUINTUPLE': rdkit.Chem.rdchem.BondType.QUINTUPLE, 'HEXTUPLE': rdkit.Chem.rdchem.BondType.HEXTUPLE, 'ONEANDAHALF': rdkit.Chem.rdchem.BondType.ONEANDAHALF, 'TWOANDAHALF': rdkit.Chem.rdchem.BondType.TWOANDAHALF, 'THREEANDAHALF': rdkit.Chem.rdchem.BondType.THREEANDAHALF, 'FOURANDAHALF': rdkit.Chem.rdchem.BondType.FOURANDAHALF, 'FIVEANDAHALF': rdkit.Chem.rdchem.BondType.FIVEANDAHALF, 'AROMATIC': rdkit.Chem.rdchem.BondType.AROMATIC, 'IONIC': rdkit.Chem.rdchem.BondType.IONIC, 'HYDROGEN': rdkit.Chem.rdchem.BondType.HYDROGEN, 'THREECENTER': rdkit.Chem.rdchem.BondType.THREECENTER, 'DATIVEONE': rdkit.Chem.rdchem.BondType.DATIVEONE, 'DATIVE': rdkit.Chem.rdchem.BondType.DATIVE, 'DATIVEL': rdkit.Chem.rdchem.BondType.DATIVEL, 'DATIVER': rdkit.Chem.rdchem.BondType.DATIVER, 'OTHER': rdkit.Chem.rdchem.BondType.OTHER, 'ZERO': rdkit.Chem.rdchem.BondType.ZERO}
    values = {0: rdkit.Chem.rdchem.BondType.UNSPECIFIED, 1: rdkit.Chem.rdchem.BondType.SINGLE, 2: rdkit.Chem.rdchem.BondType.DOUBLE, 3: rdkit.Chem.rdchem.BondType.TRIPLE, 4: rdkit.Chem.rdchem.BondType.QUADRUPLE, 5: rdkit.Chem.rdchem.BondType.QUINTUPLE, 6: rdkit.Chem.rdchem.BondType.HEXTUPLE, 7: rdkit.Chem.rdchem.BondType.ONEANDAHALF, 8: rdkit.Chem.rdchem.BondType.TWOANDAHALF, 9: rdkit.Chem.rdchem.BondType.THREEANDAHALF, 10: rdkit.Chem.rdchem.BondType.FOURANDAHALF, 11: rdkit.Chem.rdchem.BondType.FIVEANDAHALF, 12: rdkit.Chem.rdchem.BondType.AROMATIC, 13: rdkit.Chem.rdchem.BondType.IONIC, 14: rdkit.Chem.rdchem.BondType.HYDROGEN, 15: rdkit.Chem.rdchem.BondType.THREECENTER, 16: rdkit.Chem.rdchem.BondType.DATIVEONE, 17: rdkit.Chem.rdchem.BondType.DATIVE, 18: rdkit.Chem.rdchem.BondType.DATIVEL, 19: rdkit.Chem.rdchem.BondType.DATIVER, 20: rdkit.Chem.rdchem.BondType.OTHER, 21: rdkit.Chem.rdchem.BondType.ZERO}
    pass
class ChiralType(Boost.Python.enum, int):
    CHI_ALLENE = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
    CHI_OCTAHEDRAL = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
    CHI_OTHER = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
    CHI_SQUAREPLANAR = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
    CHI_TETRAHEDRAL = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
    CHI_TETRAHEDRAL_CCW = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
    CHI_TETRAHEDRAL_CW = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
    CHI_TRIGONALBIPYRAMIDAL = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
    CHI_UNSPECIFIED = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
    __slots__ = ()
    names = {'CHI_UNSPECIFIED': rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED, 'CHI_TETRAHEDRAL_CW': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, 'CHI_TETRAHEDRAL_CCW': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 'CHI_OTHER': rdkit.Chem.rdchem.ChiralType.CHI_OTHER, 'CHI_TETRAHEDRAL': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL, 'CHI_ALLENE': rdkit.Chem.rdchem.ChiralType.CHI_ALLENE, 'CHI_SQUAREPLANAR': rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR, 'CHI_TRIGONALBIPYRAMIDAL': rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL, 'CHI_OCTAHEDRAL': rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL}
    values = {0: rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED, 1: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, 2: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 3: rdkit.Chem.rdchem.ChiralType.CHI_OTHER, 4: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL, 5: rdkit.Chem.rdchem.ChiralType.CHI_ALLENE, 6: rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR, 7: rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL, 8: rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL}
    pass
class CompositeQueryType(Boost.Python.enum, int):
    COMPOSITE_AND = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
    COMPOSITE_OR = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
    COMPOSITE_XOR = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
    __slots__ = ()
    names = {'COMPOSITE_AND': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, 'COMPOSITE_OR': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR, 'COMPOSITE_XOR': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR}
    values = {0: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, 1: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR, 2: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR}
    pass
class Conformer(Boost.Python.instance):
    """
    The class to store 2D or 3D conformation of a molecule
    """
    @staticmethod
    def ClearComputedProps( arg1: Conformer) -> None: 
        """
        ClearComputedProps( arg1: Conformer) -> None
            Removes all computed properties from the conformer.
            
            

            C++ signature :
                void ClearComputedProps(RDKit::Conformer)
        """
    @staticmethod
    def ClearProp( arg1: Conformer, arg2: str) -> None: 
        """
        ClearProp( arg1: Conformer, arg2: str) -> None
            Removes a property from the conformer.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            

            C++ signature :
                void ClearProp(RDKit::Conformer,char const*)
        """
    @staticmethod
    def GetAtomPosition( arg1: Conformer, arg2: int) -> Point3D: 
        """
        GetAtomPosition( arg1: Conformer, arg2: int) -> Point3D
            Get the posistion of an atom
            

            C++ signature :
                RDGeom::Point3D GetAtomPosition(RDKit::Conformer const*,unsigned int)
        """
    @staticmethod
    def GetBoolProp( arg1: Conformer, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: Conformer, arg2: str) -> bool
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                bool GetBoolProp(RDKit::Conformer const*,char const*)
        """
    @staticmethod
    def GetDoubleProp( arg1: Conformer, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: Conformer, arg2: str) -> float
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                double GetDoubleProp(RDKit::Conformer const*,char const*)
        """
    @staticmethod
    def GetId( arg1: Conformer) -> int: 
        """
        GetId( arg1: Conformer) -> int
            Get the ID of the conformer

            C++ signature :
                unsigned int GetId(RDKit::Conformer {lvalue})
        """
    @staticmethod
    def GetIntProp( arg1: Conformer, arg2: str) -> int: 
        """
        GetIntProp( arg1: Conformer, arg2: str) -> int
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                int GetIntProp(RDKit::Conformer const*,char const*)
        """
    @staticmethod
    def GetNumAtoms( arg1: Conformer) -> int: 
        """
        GetNumAtoms( arg1: Conformer) -> int
            Get the number of atoms in the conformer
            

            C++ signature :
                unsigned int GetNumAtoms(RDKit::Conformer {lvalue})
        """
    @staticmethod
    def GetOwningMol( arg1: Conformer) -> Mol: 
        """
        GetOwningMol( arg1: Conformer) -> Mol
            Get the owning molecule
            

            C++ signature :
                RDKit::ROMol {lvalue} GetOwningMol(RDKit::Conformer {lvalue})
        """
    @staticmethod
    def GetPositions( arg1: Conformer) -> object: 
        """
        GetPositions( arg1: Conformer) -> object
            Get positions of all the atoms
            

            C++ signature :
                _object* GetPositions(RDKit::Conformer const*)
        """
    def GetProp(self, key: str, autoConvert: bool = False) -> object: 
        """
        GetProp( self: Conformer, key: str, autoConvert: bool = False) -> object
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                boost::python::api::object GetProp(RDKit::Conformer const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: Conformer, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a tuple with all property names for this conformer.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::Conformer {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: Conformer, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict
            Returns a dictionary populated with the conformer's properties.
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
              RETURNS: a dictionary
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::Conformer [,bool=False [,bool=False [,bool=True]]])
        """
    @staticmethod
    def GetUnsignedProp( arg1: Conformer, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: Conformer, arg2: str) -> int
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::Conformer const*,char const*)
        """
    @staticmethod
    def HasOwningMol( arg1: Conformer) -> bool: 
        """
        HasOwningMol( arg1: Conformer) -> bool
            Returns whether or not this instance belongs to a molecule.
            

            C++ signature :
                bool HasOwningMol(RDKit::Conformer {lvalue})
        """
    @staticmethod
    def HasProp( arg1: Conformer, arg2: str) -> int: 
        """
        HasProp( arg1: Conformer, arg2: str) -> int
            Queries a conformer to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            

            C++ signature :
                int HasProp(RDKit::Conformer,char const*)
        """
    @staticmethod
    def Is3D( arg1: Conformer) -> bool: 
        """
        Is3D( arg1: Conformer) -> bool
            returns the 3D flag of the conformer
            

            C++ signature :
                bool Is3D(RDKit::Conformer {lvalue})
        """
    @staticmethod
    def Set3D( arg1: Conformer, arg2: bool) -> None: 
        """
        Set3D( arg1: Conformer, arg2: bool) -> None
            Set the 3D flag of the conformer
            

            C++ signature :
                void Set3D(RDKit::Conformer {lvalue},bool)
        """
    @staticmethod
    @typing.overload
    def SetAtomPosition( arg1: Conformer, arg2: int, arg3: object) -> None: 
        """
        SetAtomPosition( arg1: Conformer, arg2: int, arg3: object) -> None
            Set the position of the specified atom
            

            C++ signature :
                void SetAtomPosition(RDKit::Conformer*,unsigned int,boost::python::api::object)

            C++ signature :
                void SetAtomPosition(RDKit::Conformer {lvalue},unsigned int,RDGeom::Point3D)
        """
    @staticmethod
    @typing.overload
    def SetAtomPosition( arg1: Conformer, arg2: int, arg3: Point3D) -> None: ...
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None: 
        """
        SetBoolProp( self: Conformer, key: str, val: bool, computed: bool = False) -> None
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetBoolProp(RDKit::Conformer,char const*,bool [,bool=False])
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None: 
        """
        SetDoubleProp( self: Conformer, key: str, val: float, computed: bool = False) -> None
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            

            C++ signature :
                void SetDoubleProp(RDKit::Conformer,char const*,double [,bool=False])
        """
    @staticmethod
    def SetId( arg1: Conformer, arg2: int) -> None: 
        """
        SetId( arg1: Conformer, arg2: int) -> None
            Set the ID of the conformer
            

            C++ signature :
                void SetId(RDKit::Conformer {lvalue},unsigned int)
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetIntProp( self: Conformer, key: str, val: int, computed: bool = False) -> None
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetIntProp(RDKit::Conformer,char const*,int [,bool=False])
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None: 
        """
        SetProp( self: Conformer, key: str, val: str, computed: bool = False) -> None
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetProp(RDKit::Conformer,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetUnsignedProp( self: Conformer, key: str, val: int, computed: bool = False) -> None
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetUnsignedProp(RDKit::Conformer,char const*,unsigned int [,bool=False])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,unsigned int)

            C++ signature :
                void __init__(_object*,RDKit::Conformer)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: Conformer) -> None: ...
    __instance_size__ = 40
    pass
class EditableMol(Boost.Python.instance):
    """
    an editable molecule class
    """
    @staticmethod
    def AddAtom( arg1: EditableMol, atom: Atom) -> int: 
        """
        AddAtom( arg1: EditableMol, atom: Atom) -> int
            add an atom, returns the index of the newly added atom

            C++ signature :
                int AddAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},RDKit::Atom*)
        """
    @staticmethod
    def AddBond( arg1: EditableMol, beginAtomIdx: int, endAtomIdx: int, order: BondType = BondType.UNSPECIFIED) -> int: 
        """
        AddBond( arg1: EditableMol, beginAtomIdx: int, endAtomIdx: int, order: BondType = rdkit.Chem.rdchem.BondType.UNSPECIFIED) -> int
            add a bond, returns the total number of bonds

            C++ signature :
                int AddBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,unsigned int [,RDKit::Bond::BondType=rdkit.Chem.rdchem.BondType.UNSPECIFIED])
        """
    @staticmethod
    def BeginBatchEdit( arg1: EditableMol) -> None: 
        """
        BeginBatchEdit( arg1: EditableMol) -> None
            starts batch editing

            C++ signature :
                void BeginBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
    @staticmethod
    def CommitBatchEdit( arg1: EditableMol) -> None: 
        """
        CommitBatchEdit( arg1: EditableMol) -> None
            finishes batch editing and makes the actual edits

            C++ signature :
                void CommitBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
    @staticmethod
    def GetMol( arg1: EditableMol) -> Mol: 
        """
        GetMol( arg1: EditableMol) -> Mol
            Returns a Mol (a normal molecule)

            C++ signature :
                RDKit::ROMol* GetMol(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
    @staticmethod
    def RemoveAtom( arg1: EditableMol, arg2: int) -> None: 
        """
        RemoveAtom( arg1: EditableMol, arg2: int) -> None
            Remove the specified atom from the molecule

            C++ signature :
                void RemoveAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int)
        """
    @staticmethod
    def RemoveBond( arg1: EditableMol, arg2: int, arg3: int) -> None: 
        """
        RemoveBond( arg1: EditableMol, arg2: int, arg3: int) -> None
            Remove the specified bond from the molecule

            C++ signature :
                void RemoveBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def ReplaceAtom( arg1: EditableMol, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None: 
        """
        ReplaceAtom( arg1: EditableMol, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None
            replaces the specified atom with the provided one
            If updateLabel is True, the new atom becomes the active atom
            If preserveProps is True preserve keep the existing props unless explicit set on the new atom

            C++ signature :
                void ReplaceAtom(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,RDKit::Atom* [,bool=False [,bool=False]])
        """
    @staticmethod
    def ReplaceBond( arg1: EditableMol, index: int, newBond: Bond, preserveProps: bool = False) -> None: 
        """
        ReplaceBond( arg1: EditableMol, index: int, newBond: Bond, preserveProps: bool = False) -> None
            replaces the specified bond with the provided one.
            If preserveProps is True preserve keep the existing props unless explicit set on the new bond

            C++ signature :
                void ReplaceBond(RDKit::(anonymous namespace)::EditableMol {lvalue},unsigned int,RDKit::Bond* [,bool=False])
        """
    @staticmethod
    def RollbackBatchEdit( arg1: EditableMol) -> None: 
        """
        RollbackBatchEdit( arg1: EditableMol) -> None
            cancels batch editing

            C++ signature :
                void RollbackBatchEdit(RDKit::(anonymous namespace)::EditableMol {lvalue})
        """
    @staticmethod
    def __init__( arg1: object, arg2: Mol) -> None: 
        """
        __init__( arg1: object, arg2: Mol) -> None
            Construct from a Mol

            C++ signature :
                void __init__(_object*,RDKit::ROMol)
        """
    __instance_size__ = 32
    pass
class MolBundle(Boost.Python.instance):
    """
    A class for storing groups of related molecules.
    """
    @staticmethod
    def AddMol( arg1: MolBundle, arg2: Mol) -> int: 
        """
        AddMol( arg1: MolBundle, arg2: Mol) -> int

            C++ signature :
                unsigned long AddMol(RDKit::MolBundle {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def GetMol( arg1: MolBundle, arg2: int) -> Mol: 
        """
        GetMol( arg1: MolBundle, arg2: int) -> Mol

            C++ signature :
                boost::shared_ptr<RDKit::ROMol> GetMol(RDKit::MolBundle {lvalue},unsigned long)
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: 
        """
        GetSubstructMatch( self: MolBundle, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object
            Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatch(RDKit::MolBundle,RDKit::ROMol [,bool=False [,bool=False]])

            C++ signature :
                _object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle [,bool=False [,bool=False]])

            C++ signature :
                _object* GetSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                _object* GetSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: ...
    @typing.overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: 
        """
        GetSubstructMatches( self: MolBundle, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object
            Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatches(RDKit::MolBundle,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatches(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                _object* GetSubstructMatches(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: 
        """
        HasSubstructMatch( self: MolBundle, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool
            Queries whether or not any molecule in the bundle contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - recursionPossible: (optional)
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            

            C++ signature :
                bool HasSubstructMatch(RDKit::MolBundle,RDKit::ROMol [,bool=True [,bool=False [,bool=False]]])

            C++ signature :
                bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle [,bool=True [,bool=False [,bool=False]]])

            C++ signature :
                bool HasSubstructMatch(RDKit::MolBundle,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                bool HasSubstructMatch(RDKit::MolBundle,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: ...
    @typing.overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> bool: ...
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> bool: ...
    @staticmethod
    def Size( arg1: MolBundle) -> int: 
        """
        Size( arg1: MolBundle) -> int

            C++ signature :
                unsigned long Size(RDKit::MolBundle {lvalue})
        """
    @staticmethod
    def ToBinary( arg1: MolBundle) -> object: 
        """
        ToBinary( arg1: MolBundle) -> object
            Returns a binary string representation of the MolBundle.
            

            C++ signature :
                boost::python::api::object ToBinary(RDKit::MolBundle)
        """
    @staticmethod
    def __getinitargs__( arg1: MolBundle) -> tuple: 
        """
        __getinitargs__( arg1: MolBundle) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::MolBundle)
        """
    @staticmethod
    def __getitem__( arg1: MolBundle, arg2: int) -> Mol: 
        """
        __getitem__( arg1: MolBundle, arg2: int) -> Mol

            C++ signature :
                boost::shared_ptr<RDKit::ROMol> __getitem__(RDKit::MolBundle {lvalue},unsigned long)
        """
    @staticmethod
    def __getstate__( arg1: object) -> tuple: 
        """
        __getstate__( arg1: object) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, pklString: str) -> None: ...
    @staticmethod
    def __len__( arg1: MolBundle) -> int: 
        """
        __len__( arg1: MolBundle) -> int

            C++ signature :
                unsigned long __len__(RDKit::MolBundle {lvalue})
        """
    @staticmethod
    def __setstate__( arg1: object, arg2: tuple) -> None: 
        """
        __setstate__( arg1: object, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 88
    __safe_for_unpickling__ = True
    pass
class HybridizationType(Boost.Python.enum, int):
    OTHER = rdkit.Chem.rdchem.HybridizationType.OTHER
    S = rdkit.Chem.rdchem.HybridizationType.S
    SP = rdkit.Chem.rdchem.HybridizationType.SP
    SP2 = rdkit.Chem.rdchem.HybridizationType.SP2
    SP2D = rdkit.Chem.rdchem.HybridizationType.SP2D
    SP3 = rdkit.Chem.rdchem.HybridizationType.SP3
    SP3D = rdkit.Chem.rdchem.HybridizationType.SP3D
    SP3D2 = rdkit.Chem.rdchem.HybridizationType.SP3D2
    UNSPECIFIED = rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED
    __slots__ = ()
    names = {'UNSPECIFIED': rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED, 'S': rdkit.Chem.rdchem.HybridizationType.S, 'SP': rdkit.Chem.rdchem.HybridizationType.SP, 'SP2': rdkit.Chem.rdchem.HybridizationType.SP2, 'SP3': rdkit.Chem.rdchem.HybridizationType.SP3, 'SP2D': rdkit.Chem.rdchem.HybridizationType.SP2D, 'SP3D': rdkit.Chem.rdchem.HybridizationType.SP3D, 'SP3D2': rdkit.Chem.rdchem.HybridizationType.SP3D2, 'OTHER': rdkit.Chem.rdchem.HybridizationType.OTHER}
    values = {0: rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED, 1: rdkit.Chem.rdchem.HybridizationType.S, 2: rdkit.Chem.rdchem.HybridizationType.SP, 3: rdkit.Chem.rdchem.HybridizationType.SP2, 4: rdkit.Chem.rdchem.HybridizationType.SP3, 5: rdkit.Chem.rdchem.HybridizationType.SP2D, 6: rdkit.Chem.rdchem.HybridizationType.SP3D, 7: rdkit.Chem.rdchem.HybridizationType.SP3D2, 8: rdkit.Chem.rdchem.HybridizationType.OTHER}
    pass
class KekulizeException(MolSanitizeException, ValueError, Exception, BaseException):
    pass
class Mol(Boost.Python.instance):
    """
    The Molecule class.

      In addition to the expected Atoms and Bonds, molecules contain:
        - a collection of Atom and Bond bookmarks indexed with integers
            that can be used to flag and retrieve particular Atoms or Bonds
            using the {get|set}{Atom|Bond}Bookmark() methods.

        - a set of string-valued properties. These can have arbitrary string
            labels and can be set and retrieved using the {set|get}Prop() methods
            Molecular properties can be tagged as being *computed*, in which case
              they will be automatically cleared under certain circumstances (when the
              molecule itself is modified, for example).
            Molecules also have the concept of *private* properties, which are tagged
              by beginning the property name with an underscore (_).
    """
    def AddConformer(self, conf: Conformer, assignId: bool = False) -> int: 
        """
        AddConformer( self: Mol, conf: Conformer, assignId: bool = False) -> int
            Add a conformer to the molecule and return the conformer ID

            C++ signature :
                unsigned int AddConformer(RDKit::ROMol {lvalue},RDKit::Conformer* [,bool=False])
        """
    def ClearComputedProps(self, includeRings: bool = True) -> None: 
        """
        ClearComputedProps( self: Mol, includeRings: bool = True) -> None
            Removes all computed properties from the molecule.
            
            

            C++ signature :
                void ClearComputedProps(RDKit::ROMol [,bool=True])
        """
    @staticmethod
    def ClearProp( arg1: Mol, arg2: str) -> None: 
        """
        ClearProp( arg1: Mol, arg2: str) -> None
            Removes a property from the molecule.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            

            C++ signature :
                void ClearProp(RDKit::ROMol,char const*)
        """
    def Debug(self, useStdout: bool = True) -> None: 
        """
        Debug( self: Mol, useStdout: bool = True) -> None
            Prints debugging information about the molecule.
            

            C++ signature :
                void Debug(RDKit::ROMol [,bool=True])
        """
    @staticmethod
    def GetAromaticAtoms( arg1: Mol) -> _ROQAtomSeq: 
        """
        GetAromaticAtoms( arg1: Mol) -> _ROQAtomSeq
            Returns a read-only sequence containing all of the molecule's aromatic Atoms.
            

            C++ signature :
                RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* GetAromaticAtoms(boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def GetAtomWithIdx( arg1: Mol, arg2: int) -> Atom: 
        """
        GetAtomWithIdx( arg1: Mol, arg2: int) -> Atom
            Returns a particular Atom.
            
              ARGUMENTS:
                - idx: which Atom to return
            
              NOTE: atom indices start at 0
            

            C++ signature :
                RDKit::Atom* GetAtomWithIdx(RDKit::ROMol {lvalue},unsigned int)
        """
    @staticmethod
    def GetAtomsMatchingQuery( arg1: Mol, arg2: QueryAtom) -> _ROQAtomSeq: 
        """
        GetAtomsMatchingQuery( arg1: Mol, arg2: QueryAtom) -> _ROQAtomSeq
            Returns a read-only sequence containing all of the atoms in a molecule that match the query atom.
            

            C++ signature :
                RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* GetAtomsMatchingQuery(boost::shared_ptr<RDKit::ROMol>,RDKit::QueryAtom*)
        """
    @staticmethod
    def GetBondBetweenAtoms( arg1: Mol, arg2: int, arg3: int) -> Bond: 
        """
        GetBondBetweenAtoms( arg1: Mol, arg2: int, arg3: int) -> Bond
            Returns the bond between two atoms, if there is one.
            
              ARGUMENTS:
                - idx1,idx2: the Atom indices
            
              Returns:
                The Bond between the two atoms, if such a bond exists.
                If there is no Bond between the atoms, None is returned instead.
            
              NOTE: bond indices start at 0
            

            C++ signature :
                RDKit::Bond* GetBondBetweenAtoms(RDKit::ROMol {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def GetBondWithIdx( arg1: Mol, arg2: int) -> Bond: 
        """
        GetBondWithIdx( arg1: Mol, arg2: int) -> Bond
            Returns a particular Bond.
            
              ARGUMENTS:
                - idx: which Bond to return
            
              NOTE: bond indices start at 0
            

            C++ signature :
                RDKit::Bond* GetBondWithIdx(RDKit::ROMol {lvalue},unsigned int)
        """
    @staticmethod
    def GetBoolProp( arg1: Mol, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: Mol, arg2: str) -> bool
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                bool GetBoolProp(RDKit::ROMol const*,char const*)
        """
    def GetConformer(self, id: int = -1) -> Conformer: 
        """
        GetConformer( self: Mol, id: int = -1) -> Conformer
            Get the conformer with a specified ID

            C++ signature :
                RDKit::Conformer* GetConformer(RDKit::ROMol {lvalue} [,int=-1])
        """
    @staticmethod
    def GetConformers( arg1: Mol) -> _ROConformerSeq: 
        """
        GetConformers( arg1: Mol) -> _ROConformerSeq
            Returns a read-only sequence containing all of the molecule's Conformers.

            C++ signature :
                RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor>* GetConformers(boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def GetDoubleProp( arg1: Mol, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: Mol, arg2: str) -> float
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                double GetDoubleProp(RDKit::ROMol const*,char const*)
        """
    @staticmethod
    def GetIntProp( arg1: Mol, arg2: str) -> int: 
        """
        GetIntProp( arg1: Mol, arg2: str) -> int
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                int GetIntProp(RDKit::ROMol const*,char const*)
        """
    @staticmethod
    def GetNumAtoms( arg1: Mol, onlyHeavy: int = -1, onlyExplicit: bool = True) -> int: 
        """
        GetNumAtoms( arg1: Mol, onlyHeavy: int = -1, onlyExplicit: bool = True) -> int
            Returns the number of atoms in the molecule.
            
              ARGUMENTS:
                - onlyExplicit: (optional) include only explicit atoms (atoms in the molecular graph)
                                defaults to 1.
              NOTE: the onlyHeavy argument is deprecated
            

            C++ signature :
                int GetNumAtoms(RDKit::ROMol [,int=-1 [,bool=True]])
        """
    @staticmethod
    def GetNumBonds( arg1: Mol, onlyHeavy: bool = True) -> int: 
        """
        GetNumBonds( arg1: Mol, onlyHeavy: bool = True) -> int
            Returns the number of Bonds in the molecule.
            
              ARGUMENTS:
                - onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)
                              defaults to 1.
            

            C++ signature :
                unsigned int GetNumBonds(RDKit::ROMol {lvalue} [,bool=True])
        """
    @staticmethod
    def GetNumConformers( arg1: Mol) -> int: 
        """
        GetNumConformers( arg1: Mol) -> int
            Return the number of conformations on the molecule

            C++ signature :
                unsigned int GetNumConformers(RDKit::ROMol {lvalue})
        """
    @staticmethod
    def GetNumHeavyAtoms( arg1: Mol) -> int: 
        """
        GetNumHeavyAtoms( arg1: Mol) -> int
            Returns the number of heavy atoms (atomic number >1) in the molecule.
            
            

            C++ signature :
                unsigned int GetNumHeavyAtoms(RDKit::ROMol {lvalue})
        """
    def GetProp(self, key: str, autoConvert: bool = False) -> object: 
        """
        GetProp( self: Mol, key: str, autoConvert: bool = False) -> object
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                boost::python::api::object GetProp(RDKit::ROMol const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: Mol, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a tuple with all property names for this molecule.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::ROMol {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: Mol, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict
            Returns a dictionary populated with the molecules properties.
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
              RETURNS: a dictionary
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::ROMol [,bool=False [,bool=False [,bool=True]]])
        """
    @staticmethod
    def GetRingInfo( arg1: Mol) -> RingInfo: 
        """
        GetRingInfo( arg1: Mol) -> RingInfo
            Returns the number of molecule's RingInfo object.
            
            

            C++ signature :
                RDKit::RingInfo* GetRingInfo(RDKit::ROMol {lvalue})
        """
    @staticmethod
    def GetStereoGroups( arg1: Mol) -> StereoGroup_vect: 
        """
        GetStereoGroups( arg1: Mol) -> StereoGroup_vect
            Returns a list of StereoGroups defining the relative stereochemistry of the atoms.
            

            C++ signature :
                std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > GetStereoGroups(RDKit::ROMol {lvalue})
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: 
        """
        GetSubstructMatch( self: Mol, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object
            Returns the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatch(RDKit::ROMol,RDKit::ROMol [,bool=False [,bool=False]])

            C++ signature :
                _object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,bool=False [,bool=False]])

            C++ signature :
                _object* GetSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                _object* GetSubstructMatch(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: ...
    @typing.overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: 
        """
        GetSubstructMatches( self: Mol, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object
            Returns tuples of the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatches(RDKit::ROMol,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle [,bool=True [,bool=False [,bool=False [,unsigned int=1000]]]])

            C++ signature :
                _object* GetSubstructMatches(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                _object* GetSubstructMatches(RDKit::ROMol,RDKit::MolBundle,RDKit::SubstructMatchParameters)
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters) -> object: ...
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters) -> object: ...
    @staticmethod
    def GetUnsignedProp( arg1: Mol, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: Mol, arg2: str) -> int
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::ROMol const*,char const*)
        """
    @staticmethod
    def HasProp( arg1: Mol, arg2: str) -> int: 
        """
        HasProp( arg1: Mol, arg2: str) -> int
            Queries a molecule to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            

            C++ signature :
                int HasProp(RDKit::ROMol,char const*)
        """
    @typing.overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: 
        """
        HasSubstructMatch( self: Mol, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool
            Queries whether or not the molecule contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - recursionPossible: (optional)
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            

            C++ signature :
                bool HasSubstructMatch(RDKit::ROMol,RDKit::ROMol [,bool=True [,bool=False [,bool=False]]])

            C++ signature :
                bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,bool=True [,bool=False [,bool=False]]])

            C++ signature :
                bool HasSubstructMatch(RDKit::ROMol,RDKit::ROMol,RDKit::SubstructMatchParameters)

            C++ signature :
                bool HasSubstructMatch(RDKit::ROMol,RDKit::MolBundle [,RDKit::SubstructMatchParameters=True])
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool: ...
    @typing.overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> bool: ...
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters = True) -> bool: ...
    def NeedsUpdatePropertyCache(self) -> bool: 
        """
        NeedsUpdatePropertyCache( self: Mol) -> bool
            Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
            
            

            C++ signature :
                bool NeedsUpdatePropertyCache(RDKit::ROMol {lvalue})
        """
    @staticmethod
    def RemoveAllConformers( arg1: Mol) -> None: 
        """
        RemoveAllConformers( arg1: Mol) -> None
            Remove all the conformations on the molecule

            C++ signature :
                void RemoveAllConformers(RDKit::ROMol {lvalue})
        """
    @staticmethod
    def RemoveConformer( arg1: Mol, arg2: int) -> None: 
        """
        RemoveConformer( arg1: Mol, arg2: int) -> None
            Remove the conformer with the specified ID

            C++ signature :
                void RemoveConformer(RDKit::ROMol {lvalue},unsigned int)
        """
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None: 
        """
        SetBoolProp( self: Mol, key: str, val: bool, computed: bool = False) -> None
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetBoolProp(RDKit::ROMol,char const*,bool [,bool=False])
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None: 
        """
        SetDoubleProp( self: Mol, key: str, val: float, computed: bool = False) -> None
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            

            C++ signature :
                void SetDoubleProp(RDKit::ROMol,char const*,double [,bool=False])
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetIntProp( self: Mol, key: str, val: int, computed: bool = False) -> None
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetIntProp(RDKit::ROMol,char const*,int [,bool=False])
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None: 
        """
        SetProp( self: Mol, key: str, val: str, computed: bool = False) -> None
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetProp(RDKit::ROMol,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetUnsignedProp( self: Mol, key: str, val: int, computed: bool = False) -> None
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetUnsignedProp(RDKit::ROMol,char const*,unsigned int [,bool=False])
        """
    @staticmethod
    @typing.overload
    def ToBinary( arg1: Mol) -> object: 
        """
        ToBinary( arg1: Mol) -> object
            Returns a binary string representation of the molecule.
            

            C++ signature :
                boost::python::api::object ToBinary(RDKit::ROMol)

            C++ signature :
                boost::python::api::object ToBinary(RDKit::ROMol,unsigned int)
        """
    @typing.overload
    def ToBinary(self, propertyFlags: int) -> object: ...
    def UpdatePropertyCache(self, strict: bool = True) -> None: 
        """
        UpdatePropertyCache( self: Mol, strict: bool = True) -> None
            Regenerates computed properties like implicit valence and ring information.
            
            

            C++ signature :
                void UpdatePropertyCache(RDKit::ROMol {lvalue} [,bool=True])
        """
    @staticmethod
    def __copy__( arg1: object) -> object: 
        """
        __copy__( arg1: object) -> object

            C++ signature :
                boost::python::api::object __copy__(boost::python::api::object)
        """
    @staticmethod
    def __deepcopy__( arg1: object, arg2: dict) -> object: 
        """
        __deepcopy__( arg1: object, arg2: dict) -> object

            C++ signature :
                boost::python::api::object __deepcopy__(boost::python::api::object,boost::python::dict)
        """
    @staticmethod
    def __getinitargs__( arg1: Mol) -> tuple: 
        """
        __getinitargs__( arg1: Mol) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::ROMol)
        """
    @staticmethod
    def __getstate__( arg1: object) -> tuple: 
        """
        __getstate__( arg1: object) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None
            Constructor, takes no arguments

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)

            C++ signature :
                void __init__(_object*,RDKit::ROMol [,bool=False [,int=-1]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, pklString: str) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, pklString: str, propertyFlags: int) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None: ...
    @staticmethod
    def __setstate__( arg1: object, arg2: tuple) -> None: 
        """
        __setstate__( arg1: object, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class FixedMolSizeMolBundle(MolBundle, Boost.Python.instance):
    """
    A class for storing groups of related molecules.
        Here related means that the molecules have to have the same number of atoms.
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
class AtomKekulizeException(AtomSanitizeException, MolSanitizeException, ValueError, Exception, BaseException):
    pass
class PeriodicTable(Boost.Python.instance):
    """
    A class which stores information from the Periodic Table.

      It is not possible to create a PeriodicTable object directly from Python,
      use GetPeriodicTable() to get the global table.

      The PeriodicTable object can be queried for a variety of properties:

        - GetAtomicWeight

        - GetAtomicNumber

        - GetElementSymbol

        - GetElementName

        - GetRvdw (van der Waals radius)

        - GetRCovalent (covalent radius)

        - GetDefaultValence

        - GetValenceList

        - GetNOuterElecs (number of valence electrons)

        - GetMostCommonIsotope

        - GetMostCommonIsotopeMass

        - GetRb0

        - GetAbundanceForIsotope

        - GetMassForIsotope

      When it makes sense, these can be queried using either an atomic number (integer)
      or an atomic symbol (string)
    """
    @staticmethod
    @typing.overload
    def GetAbundanceForIsotope( arg1: PeriodicTable, arg2: int, arg3: int) -> float: 
        """
        GetAbundanceForIsotope( arg1: PeriodicTable, arg2: int, arg3: int) -> float

            C++ signature :
                double GetAbundanceForIsotope(RDKit::PeriodicTable {lvalue},unsigned int,unsigned int)

            C++ signature :
                double GetAbundanceForIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
        """
    @staticmethod
    @typing.overload
    def GetAbundanceForIsotope( arg1: PeriodicTable, arg2: str, arg3: int) -> float: ...
    @staticmethod
    def GetAtomicNumber( arg1: PeriodicTable, arg2: str) -> int: 
        """
        GetAtomicNumber( arg1: PeriodicTable, arg2: str) -> int

            C++ signature :
                int GetAtomicNumber(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetAtomicWeight( arg1: PeriodicTable, arg2: int) -> float: 
        """
        GetAtomicWeight( arg1: PeriodicTable, arg2: int) -> float

            C++ signature :
                double GetAtomicWeight(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                double GetAtomicWeight(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetAtomicWeight( arg1: PeriodicTable, arg2: str) -> float: ...
    @staticmethod
    @typing.overload
    def GetDefaultValence( arg1: PeriodicTable, arg2: int) -> int: 
        """
        GetDefaultValence( arg1: PeriodicTable, arg2: int) -> int

            C++ signature :
                int GetDefaultValence(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                int GetDefaultValence(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetDefaultValence( arg1: PeriodicTable, arg2: str) -> int: ...
    @staticmethod
    def GetElementName( arg1: PeriodicTable, arg2: int) -> str: 
        """
        GetElementName( arg1: PeriodicTable, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetElementName(RDKit::PeriodicTable {lvalue},unsigned int)
        """
    @staticmethod
    def GetElementSymbol( arg1: PeriodicTable, arg2: int) -> str: 
        """
        GetElementSymbol( arg1: PeriodicTable, arg2: int) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetElementSymbol(RDKit::PeriodicTable {lvalue},unsigned int)
        """
    @staticmethod
    @typing.overload
    def GetMassForIsotope( arg1: PeriodicTable, arg2: int, arg3: int) -> float: 
        """
        GetMassForIsotope( arg1: PeriodicTable, arg2: int, arg3: int) -> float

            C++ signature :
                double GetMassForIsotope(RDKit::PeriodicTable {lvalue},unsigned int,unsigned int)

            C++ signature :
                double GetMassForIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
        """
    @staticmethod
    @typing.overload
    def GetMassForIsotope( arg1: PeriodicTable, arg2: str, arg3: int) -> float: ...
    @staticmethod
    @typing.overload
    def GetMostCommonIsotope( arg1: PeriodicTable, arg2: int) -> int: 
        """
        GetMostCommonIsotope( arg1: PeriodicTable, arg2: int) -> int

            C++ signature :
                int GetMostCommonIsotope(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                int GetMostCommonIsotope(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetMostCommonIsotope( arg1: PeriodicTable, arg2: str) -> int: ...
    @staticmethod
    @typing.overload
    def GetMostCommonIsotopeMass( arg1: PeriodicTable, arg2: int) -> float: 
        """
        GetMostCommonIsotopeMass( arg1: PeriodicTable, arg2: int) -> float

            C++ signature :
                double GetMostCommonIsotopeMass(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                double GetMostCommonIsotopeMass(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetMostCommonIsotopeMass( arg1: PeriodicTable, arg2: str) -> float: ...
    @staticmethod
    @typing.overload
    def GetNOuterElecs( arg1: PeriodicTable, arg2: int) -> int: 
        """
        GetNOuterElecs( arg1: PeriodicTable, arg2: int) -> int

            C++ signature :
                int GetNOuterElecs(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                int GetNOuterElecs(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetNOuterElecs( arg1: PeriodicTable, arg2: str) -> int: ...
    @staticmethod
    @typing.overload
    def GetRb0( arg1: PeriodicTable, arg2: int) -> float: 
        """
        GetRb0( arg1: PeriodicTable, arg2: int) -> float

            C++ signature :
                double GetRb0(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                double GetRb0(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetRb0( arg1: PeriodicTable, arg2: str) -> float: ...
    @staticmethod
    @typing.overload
    def GetRcovalent( arg1: PeriodicTable, arg2: int) -> float: 
        """
        GetRcovalent( arg1: PeriodicTable, arg2: int) -> float

            C++ signature :
                double GetRcovalent(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                double GetRcovalent(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetRcovalent( arg1: PeriodicTable, arg2: str) -> float: ...
    @staticmethod
    @typing.overload
    def GetRvdw( arg1: PeriodicTable, arg2: int) -> float: 
        """
        GetRvdw( arg1: PeriodicTable, arg2: int) -> float

            C++ signature :
                double GetRvdw(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                double GetRvdw(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetRvdw( arg1: PeriodicTable, arg2: str) -> float: ...
    @staticmethod
    @typing.overload
    def GetValenceList( arg1: PeriodicTable, arg2: int) -> _vecti: 
        """
        GetValenceList( arg1: PeriodicTable, arg2: int) -> _vecti

            C++ signature :
                std::vector<int, std::allocator<int> > GetValenceList(RDKit::PeriodicTable {lvalue},unsigned int)

            C++ signature :
                std::vector<int, std::allocator<int> > GetValenceList(RDKit::PeriodicTable {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def GetValenceList( arg1: PeriodicTable, arg2: str) -> _vecti: ...
    pass
class PropertyPickleOptions(Boost.Python.enum, int):
    AllProps = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
    AtomProps = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
    BondProps = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
    ComputedProps = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
    CoordsAsDouble = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
    MolProps = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
    NoConformers = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
    NoProps = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
    PrivateProps = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
    QueryAtomData = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
    __slots__ = ()
    names = {'NoProps': rdkit.Chem.rdchem.PropertyPickleOptions.NoProps, 'MolProps': rdkit.Chem.rdchem.PropertyPickleOptions.MolProps, 'AtomProps': rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps, 'BondProps': rdkit.Chem.rdchem.PropertyPickleOptions.BondProps, 'QueryAtomData': rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData, 'PrivateProps': rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps, 'ComputedProps': rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps, 'AllProps': rdkit.Chem.rdchem.PropertyPickleOptions.AllProps, 'CoordsAsDouble': rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble, 'NoConformers': rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers}
    values = {0: rdkit.Chem.rdchem.PropertyPickleOptions.NoProps, 1: rdkit.Chem.rdchem.PropertyPickleOptions.MolProps, 2: rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData, 4: rdkit.Chem.rdchem.PropertyPickleOptions.BondProps, 16: rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps, 32: rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps, 65535: rdkit.Chem.rdchem.PropertyPickleOptions.AllProps, 65536: rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble, 131072: rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers}
    pass
class QueryAtom(Atom, Boost.Python.instance):
    """
    The class to store QueryAtoms.
    These cannot currently be constructed directly from Python
    """
    def ExpandQuery(self, other: QueryAtom, how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None: 
        """
        ExpandQuery( self: QueryAtom, other: QueryAtom, how: CompositeQueryType = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None
            combines the query from other with ours

            C++ signature :
                void ExpandQuery(RDKit::QueryAtom*,RDKit::QueryAtom const* [,Queries::CompositeQueryType=rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND [,bool=True]])
        """
    def SetQuery(self, other: QueryAtom) -> None: 
        """
        SetQuery( self: QueryAtom, other: QueryAtom) -> None
            Replace our query with a copy of the other query

            C++ signature :
                void SetQuery(RDKit::QueryAtom*,RDKit::QueryAtom const*)
        """
    pass
class QueryBond(Bond, Boost.Python.instance):
    """
    The class to store QueryBonds.
    These cannot currently be constructed directly from Python
    """
    def ExpandQuery(self, other: QueryBond, how: CompositeQueryType = CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None: 
        """
        ExpandQuery( self: QueryBond, other: QueryBond, how: CompositeQueryType = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, maintainOrder: bool = True) -> None
            combines the query from other with ours

            C++ signature :
                void ExpandQuery(RDKit::QueryBond*,RDKit::QueryBond const* [,Queries::CompositeQueryType=rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND [,bool=True]])
        """
    def SetQuery(self, other: QueryBond) -> None: 
        """
        SetQuery( self: QueryBond, other: QueryBond) -> None
            Replace our query with a copy of the other query

            C++ signature :
                void SetQuery(RDKit::QueryBond*,RDKit::QueryBond const*)
        """
    pass
class RWMol(Mol, Boost.Python.instance):
    """
    The RW molecule class (read/write)

      This class is a more-performant version of the EditableMolecule class in that
      it is a 'live' molecule and shares the interface from the Mol class.
      All changes are performed without the need to create a copy of the
      molecule using GetMol() (this is still available, however).
      
      n.b. Eventually this class may become a direct replacement for EditableMol
    """
    @staticmethod
    def AddAtom( arg1: RWMol, atom: Atom) -> int: 
        """
        AddAtom( arg1: RWMol, atom: Atom) -> int
            add an atom, returns the index of the newly added atom

            C++ signature :
                int AddAtom(RDKit::ReadWriteMol {lvalue},RDKit::Atom*)
        """
    @staticmethod
    def AddBond( arg1: RWMol, beginAtomIdx: int, endAtomIdx: int, order: BondType = BondType.UNSPECIFIED) -> int: 
        """
        AddBond( arg1: RWMol, beginAtomIdx: int, endAtomIdx: int, order: BondType = rdkit.Chem.rdchem.BondType.UNSPECIFIED) -> int
            add a bond, returns the new number of bonds

            C++ signature :
                int AddBond(RDKit::ReadWriteMol {lvalue},unsigned int,unsigned int [,RDKit::Bond::BondType=rdkit.Chem.rdchem.BondType.UNSPECIFIED])
        """
    @staticmethod
    def BeginBatchEdit( arg1: RWMol) -> None: 
        """
        BeginBatchEdit( arg1: RWMol) -> None
            starts batch editing

            C++ signature :
                void BeginBatchEdit(RDKit::ReadWriteMol {lvalue})
        """
    @staticmethod
    def CommitBatchEdit( arg1: RWMol) -> None: 
        """
        CommitBatchEdit( arg1: RWMol) -> None
            finishes batch editing and makes the actual changes

            C++ signature :
                void CommitBatchEdit(RDKit::ReadWriteMol {lvalue})
        """
    @staticmethod
    def GetMol( arg1: RWMol) -> Mol: 
        """
        GetMol( arg1: RWMol) -> Mol
            Returns a Mol (a normal molecule)

            C++ signature :
                RDKit::ROMol* GetMol(RDKit::ReadWriteMol {lvalue})
        """
    @staticmethod
    def InsertMol( arg1: RWMol, mol: Mol) -> None: 
        """
        InsertMol( arg1: RWMol, mol: Mol) -> None
            Insert (add) the given molecule into this one

            C++ signature :
                void InsertMol(RDKit::ReadWriteMol {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def RemoveAtom( arg1: RWMol, arg2: int) -> None: 
        """
        RemoveAtom( arg1: RWMol, arg2: int) -> None
            Remove the specified atom from the molecule

            C++ signature :
                void RemoveAtom(RDKit::ReadWriteMol {lvalue},unsigned int)
        """
    @staticmethod
    def RemoveBond( arg1: RWMol, arg2: int, arg3: int) -> None: 
        """
        RemoveBond( arg1: RWMol, arg2: int, arg3: int) -> None
            Remove the specified bond from the molecule

            C++ signature :
                void RemoveBond(RDKit::ReadWriteMol {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def ReplaceAtom( arg1: RWMol, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None: 
        """
        ReplaceAtom( arg1: RWMol, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None
            replaces the specified atom with the provided one
            If updateLabel is True, the new atom becomes the active atom
            If preserveProps is True preserve keep the existing props unless explicit set on the new atom

            C++ signature :
                void ReplaceAtom(RDKit::ReadWriteMol {lvalue},unsigned int,RDKit::Atom* [,bool=False [,bool=False]])
        """
    @staticmethod
    def ReplaceBond( arg1: RWMol, index: int, newBond: Bond, preserveProps: bool = False, keepSGroups: bool = True) -> None: 
        """
        ReplaceBond( arg1: RWMol, index: int, newBond: Bond, preserveProps: bool = False, keepSGroups: bool = True) -> None
            replaces the specified bond with the provided one.
            If preserveProps is True preserve keep the existing props unless explicit set on the new bond. If keepSGroups is False, allSubstance Groups referencing the bond will be dropped.

            C++ signature :
                void ReplaceBond(RDKit::ReadWriteMol {lvalue},unsigned int,RDKit::Bond* [,bool=False [,bool=True]])
        """
    @staticmethod
    def RollbackBatchEdit( arg1: RWMol) -> None: 
        """
        RollbackBatchEdit( arg1: RWMol) -> None
            cancels batch editing

            C++ signature :
                void RollbackBatchEdit(RDKit::ReadWriteMol {lvalue})
        """
    @staticmethod
    def SetStereoGroups( arg1: RWMol, stereo_groups: list) -> None: 
        """
        SetStereoGroups( arg1: RWMol, stereo_groups: list) -> None
            Set the stereo groups

            C++ signature :
                void SetStereoGroups(RDKit::ReadWriteMol {lvalue},boost::python::list {lvalue})
        """
    @staticmethod
    def __copy__( arg1: object) -> object: 
        """
        __copy__( arg1: object) -> object

            C++ signature :
                boost::python::api::object __copy__(boost::python::api::object)
        """
    @staticmethod
    def __deepcopy__( arg1: object, arg2: dict) -> object: 
        """
        __deepcopy__( arg1: object, arg2: dict) -> object

            C++ signature :
                boost::python::api::object __deepcopy__(boost::python::api::object,boost::python::dict)
        """
    @staticmethod
    def __enter__( arg1: RWMol) -> RWMol: 
        """
        __enter__( arg1: RWMol) -> RWMol

            C++ signature :
                RDKit::ReadWriteMol* __enter__(RDKit::ReadWriteMol {lvalue})
        """
    @staticmethod
    def __exit__( arg1: RWMol, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: RWMol, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::ReadWriteMol {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def __getinitargs__( arg1: Mol) -> tuple: 
        """
        __getinitargs__( arg1: Mol) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::ROMol)
        """
    @staticmethod
    def __getstate__( arg1: object) -> tuple: 
        """
        __getstate__( arg1: object) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: Mol) -> None: 
        """
        __init__( arg1: object, arg2: Mol) -> None
            Construct from a Mol

            C++ signature :
                void __init__(_object*,RDKit::ROMol)

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)

            C++ signature :
                void __init__(_object*,RDKit::ROMol [,bool=False [,int=-1]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, pklString: str) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, pklString: str, propertyFlags: int) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None: ...
    @staticmethod
    def __setstate__( arg1: object, arg2: tuple) -> None: 
        """
        __setstate__( arg1: object, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 336
    __safe_for_unpickling__ = True
    pass
class ResonanceFlags(Boost.Python.enum, int):
    ALLOW_CHARGE_SEPARATION = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
    ALLOW_INCOMPLETE_OCTETS = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
    KEKULE_ALL = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
    UNCONSTRAINED_ANIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
    UNCONSTRAINED_CATIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
    __slots__ = ()
    names = {'ALLOW_INCOMPLETE_OCTETS': rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS, 'ALLOW_CHARGE_SEPARATION': rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION, 'KEKULE_ALL': rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL, 'UNCONSTRAINED_CATIONS': rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS, 'UNCONSTRAINED_ANIONS': rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS}
    values = {1: rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS, 2: rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION, 4: rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL, 8: rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS, 16: rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS}
    pass
class ResonanceMolSupplier(Boost.Python.instance):
    """
    A class which supplies resonance structures (as mols) from a mol.

      Usage examples:

        1) Lazy evaluation: the resonance structures are not constructed
           until we ask for them:

           >>> suppl = ResonanceMolSupplier(mol)
           >>> for resMol in suppl:
           ...    resMol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = ResonanceMolSupplier(mol)
           >>> resMol1 = next(suppl)
           >>> resMol2 = next(suppl)
           >>> suppl.reset()
           >>> resMol3 = next(suppl)
           # resMol3 and resMol1 are the same: 
           >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)

        3) Random Access:

           >>> suppl = ResonanceMolSupplier(mol)
           >>> resMol1 = suppl[0] 
           >>> resMol2 = suppl[1] 

           NOTE: this will generate an IndexError if the supplier doesn't have that many
           molecules.

        4) Random Access 2: looping over all resonance structures
           >>> suppl = ResonanceMolSupplier(mol)
           >>> nResMols = len(suppl)
           >>> for i in range(nResMols):
           ...   suppl[i].GetNumAtoms()
    """
    @staticmethod
    def Enumerate( arg1: ResonanceMolSupplier) -> None: 
        """
        Enumerate( arg1: ResonanceMolSupplier) -> None
            Ask ResonanceMolSupplier to enumerate resonance structures(automatically done as soon as any attempt to access them is made).
            

            C++ signature :
                void Enumerate(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def GetAtomConjGrpIdx( arg1: ResonanceMolSupplier, arg2: int) -> int: 
        """
        GetAtomConjGrpIdx( arg1: ResonanceMolSupplier, arg2: int) -> int
            Given an atom index, it returns the index of the conjugated groupthe atom belongs to, or -1 if it is not conjugated.
            

            C++ signature :
                unsigned int GetAtomConjGrpIdx(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
    @staticmethod
    def GetBondConjGrpIdx( arg1: ResonanceMolSupplier, arg2: int) -> int: 
        """
        GetBondConjGrpIdx( arg1: ResonanceMolSupplier, arg2: int) -> int
            Given a bond index, it returns the index of the conjugated groupthe bond belongs to, or -1 if it is not conjugated.
            

            C++ signature :
                unsigned int GetBondConjGrpIdx(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
    @staticmethod
    def GetIsEnumerated( arg1: ResonanceMolSupplier) -> bool: 
        """
        GetIsEnumerated( arg1: ResonanceMolSupplier) -> bool
            Returns true if resonance structure enumeration has already happened.
            

            C++ signature :
                bool GetIsEnumerated(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def GetNumConjGrps( arg1: ResonanceMolSupplier) -> int: 
        """
        GetNumConjGrps( arg1: ResonanceMolSupplier) -> int
            Returns the number of individual conjugated groups in the molecule.
            

            C++ signature :
                unsigned int GetNumConjGrps(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def GetProgressCallback( arg1: ResonanceMolSupplier) -> object: 
        """
        GetProgressCallback( arg1: ResonanceMolSupplier) -> object
            Get the ResonanceMolSupplierCallback subclass instance,
            or None if none was set.
            

            C++ signature :
                boost::python::api::object GetProgressCallback(RDKit::ResonanceMolSupplier)
        """
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object: 
        """
        GetSubstructMatch( self: ResonanceMolSupplier, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> object
            Returns the indices of the molecule's atoms that match a substructure query,
            taking into account all resonance structures in ResonanceMolSupplier.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatch(RDKit::ResonanceMolSupplier {lvalue},RDKit::ROMol [,bool=False [,bool=False]])
        """
    def GetSubstructMatches(self, query: Mol, uniquify: bool = False, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000, numThreads: int = 1) -> object: 
        """
        GetSubstructMatches( self: ResonanceMolSupplier, query: Mol, uniquify: bool = False, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000, numThreads: int = 1) -> object
            Returns tuples of the indices of the molecule's atoms that match a substructure query,
            taking into account all resonance structures in ResonanceMolSupplier.
            
              ARGUMENTS:
                - query: a Molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
                - numThreads: The number of threads to be used (defaults to 1; 0 selects the
                              number of concurrent threads supported by the hardware; negative
                              values are added to the number of concurrent threads supported
                              by the hardware).
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            

            C++ signature :
                _object* GetSubstructMatches(RDKit::ResonanceMolSupplier {lvalue},RDKit::ROMol [,bool=False [,bool=False [,bool=False [,unsigned int=1000 [,int=1]]]]])
        """
    @staticmethod
    def SetNumThreads( arg1: ResonanceMolSupplier, arg2: int) -> None: 
        """
        SetNumThreads( arg1: ResonanceMolSupplier, arg2: int) -> None
            Sets the number of threads to be used to enumerate resonance
            structures (defaults to 1; 0 selects the number of concurrent
            threads supported by the hardware; negative values are added
            to the number of concurrent threads supported by the hardware).
            

            C++ signature :
                void SetNumThreads(RDKit::ResonanceMolSupplier {lvalue},unsigned int)
        """
    @staticmethod
    def SetProgressCallback( arg1: ResonanceMolSupplier, arg2: object) -> None: 
        """
        SetProgressCallback( arg1: ResonanceMolSupplier, arg2: object) -> None
            Pass an instance of a class derived from
            ResonanceMolSupplierCallback, which must implement the
            __call__() method.
            

            C++ signature :
                void SetProgressCallback(RDKit::ResonanceMolSupplier {lvalue},_object*)
        """
    @staticmethod
    def WasCanceled( arg1: ResonanceMolSupplier) -> bool: 
        """
        WasCanceled( arg1: ResonanceMolSupplier) -> bool
            Returns True if the resonance structure generation was canceled.
            

            C++ signature :
                bool WasCanceled(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def __getitem__( arg1: ResonanceMolSupplier, arg2: int) -> Mol: 
        """
        __getitem__( arg1: ResonanceMolSupplier, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__(RDKit::ResonanceMolSupplier*,int)
        """
    @staticmethod
    def __init__( arg1: object, mol: Mol, flags: int = 0, maxStructs: int = 1000) -> None: 
        """
        __init__( arg1: object, mol: Mol, flags: int = 0, maxStructs: int = 1000) -> None

            C++ signature :
                void __init__(_object*,RDKit::ROMol {lvalue} [,unsigned int=0 [,unsigned int=1000]])
        """
    @staticmethod
    def __iter__( arg1: ResonanceMolSupplier) -> ResonanceMolSupplier: 
        """
        __iter__( arg1: ResonanceMolSupplier) -> ResonanceMolSupplier

            C++ signature :
                RDKit::ResonanceMolSupplier* __iter__(RDKit::ResonanceMolSupplier*)
        """
    @staticmethod
    def __len__( arg1: ResonanceMolSupplier) -> int: 
        """
        __len__( arg1: ResonanceMolSupplier) -> int

            C++ signature :
                unsigned int __len__(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def __next__( arg1: ResonanceMolSupplier) -> Mol: 
        """
        __next__( arg1: ResonanceMolSupplier) -> Mol
            Returns the next resonance structure in the supplier. Raises _StopIteration_ on end.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::ResonanceMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: ResonanceMolSupplier) -> bool: 
        """
        atEnd( arg1: ResonanceMolSupplier) -> bool
            Returns whether or not we have hit the end of the resonance structure supplier.
            

            C++ signature :
                bool atEnd(RDKit::ResonanceMolSupplier {lvalue})
        """
    @staticmethod
    def reset( arg1: ResonanceMolSupplier) -> None: 
        """
        reset( arg1: ResonanceMolSupplier) -> None
            Resets our position in the resonance structure supplier to the beginning.
            

            C++ signature :
                void reset(RDKit::ResonanceMolSupplier {lvalue})
        """
    __instance_size__ = 168
    pass
class ResonanceMolSupplierCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
        implement the __call__() method.
        The __call__() method is called at each iteration of the
        algorithm, and provides a mechanism to monitor or stop
        its progress.

        To have your callback called, pass an instance of your
        derived class to ResonanceMolSupplier.SetProgressCallback()
    """
    @staticmethod
    def GetMaxStructures( arg1: ResonanceMolSupplierCallback) -> int: 
        """
        GetMaxStructures( arg1: ResonanceMolSupplierCallback) -> int
            Get the number of conjugated groups this molecule has.
            

            C++ signature :
                unsigned long GetMaxStructures(RDKit::PyResonanceMolSupplierCallback {lvalue})
        """
    @staticmethod
    def GetNumConjGrps( arg1: ResonanceMolSupplierCallback) -> int: 
        """
        GetNumConjGrps( arg1: ResonanceMolSupplierCallback) -> int
            Returns the number of individual conjugated groups in the molecule.
            

            C++ signature :
                unsigned int GetNumConjGrps(RDKit::PyResonanceMolSupplierCallback {lvalue})
        """
    @staticmethod
    def GetNumDiverseStructures( arg1: ResonanceMolSupplierCallback, arg2: int) -> int: 
        """
        GetNumDiverseStructures( arg1: ResonanceMolSupplierCallback, arg2: int) -> int
            Get the number of non-degenrate resonance structures generated so far for the passed conjugated group index.
            

            C++ signature :
                unsigned long GetNumDiverseStructures(RDKit::PyResonanceMolSupplierCallback {lvalue},unsigned int)
        """
    @staticmethod
    def GetNumStructures( arg1: ResonanceMolSupplierCallback, arg2: int) -> int: 
        """
        GetNumStructures( arg1: ResonanceMolSupplierCallback, arg2: int) -> int
            Get the number of resonance structures generated so far for the passed conjugated group index.
            

            C++ signature :
                unsigned long GetNumStructures(RDKit::PyResonanceMolSupplierCallback {lvalue},unsigned int)
        """
    @staticmethod
    @typing.overload
    def __call__( arg1: ResonanceMolSupplierCallback) -> bool: 
        """
        __call__( arg1: ResonanceMolSupplierCallback) -> bool
            This must be implemented in the derived class. Return True if the resonance structure generation should continue; False if the resonance structure generation should stop.
            

            C++ signature :
                bool __call__(RDKit::PyResonanceMolSupplierCallback {lvalue})

            C++ signature :
                void __call__(RDKit::PyResonanceMolSupplierCallback {lvalue})
        """
    @staticmethod
    @typing.overload
    def __call__( arg1: ResonanceMolSupplierCallback) -> None: ...
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 88
    pass
class RingInfo(Boost.Python.instance):
    """
    contains information about a molecule's rings
    """
    def AddRing(self, atomIds: object, bondIds: object) -> None: 
        """
        AddRing( self: RingInfo, atomIds: object, bondIds: object) -> None
            Adds a ring to the set. Be very careful with this operation.

            C++ signature :
                void AddRing(RDKit::RingInfo*,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def AreAtomsInSameRing( arg1: RingInfo, arg2: int, arg3: int) -> bool: 
        """
        AreAtomsInSameRing( arg1: RingInfo, arg2: int, arg3: int) -> bool

            C++ signature :
                bool AreAtomsInSameRing(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def AreAtomsInSameRingOfSize( arg1: RingInfo, arg2: int, arg3: int, arg4: int) -> bool: 
        """
        AreAtomsInSameRingOfSize( arg1: RingInfo, arg2: int, arg3: int, arg4: int) -> bool

            C++ signature :
                bool AreAtomsInSameRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int,unsigned int)
        """
    @staticmethod
    def AreBondsInSameRing( arg1: RingInfo, arg2: int, arg3: int) -> bool: 
        """
        AreBondsInSameRing( arg1: RingInfo, arg2: int, arg3: int) -> bool

            C++ signature :
                bool AreBondsInSameRing(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def AreBondsInSameRingOfSize( arg1: RingInfo, arg2: int, arg3: int, arg4: int) -> bool: 
        """
        AreBondsInSameRingOfSize( arg1: RingInfo, arg2: int, arg3: int, arg4: int) -> bool

            C++ signature :
                bool AreBondsInSameRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int,unsigned int)
        """
    @staticmethod
    def AreRingFamiliesInitialized( arg1: RingInfo) -> bool: 
        """
        AreRingFamiliesInitialized( arg1: RingInfo) -> bool

            C++ signature :
                bool AreRingFamiliesInitialized(RDKit::RingInfo {lvalue})
        """
    @staticmethod
    def AreRingsFused( arg1: RingInfo, arg2: int, arg3: int) -> bool: 
        """
        AreRingsFused( arg1: RingInfo, arg2: int, arg3: int) -> bool

            C++ signature :
                bool AreRingsFused(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def AtomMembers( arg1: RingInfo, arg2: int) -> object: 
        """
        AtomMembers( arg1: RingInfo, arg2: int) -> object

            C++ signature :
                boost::python::api::object AtomMembers(RDKit::RingInfo const*,unsigned int)
        """
    @staticmethod
    def AtomRingFamilies( arg1: RingInfo) -> object: 
        """
        AtomRingFamilies( arg1: RingInfo) -> object

            C++ signature :
                boost::python::api::object AtomRingFamilies(RDKit::RingInfo const*)
        """
    @staticmethod
    def AtomRingSizes( arg1: RingInfo, arg2: int) -> object: 
        """
        AtomRingSizes( arg1: RingInfo, arg2: int) -> object

            C++ signature :
                boost::python::api::object AtomRingSizes(RDKit::RingInfo const*,unsigned int)
        """
    @staticmethod
    def AtomRings( arg1: RingInfo) -> object: 
        """
        AtomRings( arg1: RingInfo) -> object

            C++ signature :
                boost::python::api::object AtomRings(RDKit::RingInfo const*)
        """
    @staticmethod
    def BondMembers( arg1: RingInfo, arg2: int) -> object: 
        """
        BondMembers( arg1: RingInfo, arg2: int) -> object

            C++ signature :
                boost::python::api::object BondMembers(RDKit::RingInfo const*,unsigned int)
        """
    @staticmethod
    def BondRingFamilies( arg1: RingInfo) -> object: 
        """
        BondRingFamilies( arg1: RingInfo) -> object

            C++ signature :
                boost::python::api::object BondRingFamilies(RDKit::RingInfo const*)
        """
    @staticmethod
    def BondRingSizes( arg1: RingInfo, arg2: int) -> object: 
        """
        BondRingSizes( arg1: RingInfo, arg2: int) -> object

            C++ signature :
                boost::python::api::object BondRingSizes(RDKit::RingInfo const*,unsigned int)
        """
    @staticmethod
    def BondRings( arg1: RingInfo) -> object: 
        """
        BondRings( arg1: RingInfo) -> object

            C++ signature :
                boost::python::api::object BondRings(RDKit::RingInfo const*)
        """
    @staticmethod
    def IsAtomInRingOfSize( arg1: RingInfo, arg2: int, arg3: int) -> bool: 
        """
        IsAtomInRingOfSize( arg1: RingInfo, arg2: int, arg3: int) -> bool

            C++ signature :
                bool IsAtomInRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def IsBondInRingOfSize( arg1: RingInfo, arg2: int, arg3: int) -> bool: 
        """
        IsBondInRingOfSize( arg1: RingInfo, arg2: int, arg3: int) -> bool

            C++ signature :
                bool IsBondInRingOfSize(RDKit::RingInfo {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def IsRingFused( arg1: RingInfo, arg2: int) -> bool: 
        """
        IsRingFused( arg1: RingInfo, arg2: int) -> bool

            C++ signature :
                bool IsRingFused(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def MinAtomRingSize( arg1: RingInfo, arg2: int) -> int: 
        """
        MinAtomRingSize( arg1: RingInfo, arg2: int) -> int

            C++ signature :
                unsigned int MinAtomRingSize(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def MinBondRingSize( arg1: RingInfo, arg2: int) -> int: 
        """
        MinBondRingSize( arg1: RingInfo, arg2: int) -> int

            C++ signature :
                unsigned int MinBondRingSize(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def NumAtomRings( arg1: RingInfo, arg2: int) -> int: 
        """
        NumAtomRings( arg1: RingInfo, arg2: int) -> int

            C++ signature :
                unsigned int NumAtomRings(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def NumBondRings( arg1: RingInfo, arg2: int) -> int: 
        """
        NumBondRings( arg1: RingInfo, arg2: int) -> int

            C++ signature :
                unsigned int NumBondRings(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def NumFusedBonds( arg1: RingInfo, arg2: int) -> int: 
        """
        NumFusedBonds( arg1: RingInfo, arg2: int) -> int

            C++ signature :
                unsigned int NumFusedBonds(RDKit::RingInfo {lvalue},unsigned int)
        """
    @staticmethod
    def NumRelevantCycles( arg1: RingInfo) -> int: 
        """
        NumRelevantCycles( arg1: RingInfo) -> int

            C++ signature :
                unsigned int NumRelevantCycles(RDKit::RingInfo {lvalue})
        """
    @staticmethod
    def NumRingFamilies( arg1: RingInfo) -> int: 
        """
        NumRingFamilies( arg1: RingInfo) -> int

            C++ signature :
                unsigned int NumRingFamilies(RDKit::RingInfo {lvalue})
        """
    @staticmethod
    def NumRings( arg1: RingInfo) -> int: 
        """
        NumRings( arg1: RingInfo) -> int

            C++ signature :
                unsigned int NumRings(RDKit::RingInfo {lvalue})
        """
    pass
class StereoDescriptor(Boost.Python.enum, int):
    Bond_Cis = rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis
    Bond_Trans = rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans
    NoValue = rdkit.Chem.rdchem.StereoDescriptor.NoValue
    Tet_CCW = rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW
    Tet_CW = rdkit.Chem.rdchem.StereoDescriptor.Tet_CW
    __slots__ = ()
    names = {'NoValue': rdkit.Chem.rdchem.StereoDescriptor.NoValue, 'Tet_CW': rdkit.Chem.rdchem.StereoDescriptor.Tet_CW, 'Tet_CCW': rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW, 'Bond_Cis': rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis, 'Bond_Trans': rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans}
    values = {0: rdkit.Chem.rdchem.StereoDescriptor.NoValue, 1: rdkit.Chem.rdchem.StereoDescriptor.Tet_CW, 2: rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW, 3: rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis, 4: rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans}
    pass
class StereoGroup(Boost.Python.instance):
    """
    A collection of atoms with a defined stereochemical relationship.

    Used to help represent a sample with unknown stereochemistry, or that is a mix
    of diastereomers.
    """
    @staticmethod
    def GetAtoms( arg1: StereoGroup) -> object: 
        """
        GetAtoms( arg1: StereoGroup) -> object
            access the atoms in the StereoGroup.
            

            C++ signature :
                boost::python::api::object GetAtoms(RDKit::StereoGroup {lvalue})
        """
    @staticmethod
    def GetGroupType( arg1: StereoGroup) -> StereoGroupType: 
        """
        GetGroupType( arg1: StereoGroup) -> StereoGroupType
            Returns the StereoGroupType.
            

            C++ signature :
                RDKit::StereoGroupType GetGroupType(RDKit::StereoGroup {lvalue})
        """
    @staticmethod
    def GetReadId( arg1: StereoGroup) -> int: 
        """
        GetReadId( arg1: StereoGroup) -> int
            return the StereoGroup's original ID.
            Note that the ID only makes sense for AND/OR groups.
            

            C++ signature :
                unsigned int GetReadId(RDKit::StereoGroup {lvalue})
        """
    @staticmethod
    def GetWriteId( arg1: StereoGroup) -> int: 
        """
        GetWriteId( arg1: StereoGroup) -> int
            return the StereoGroup's ID that will be exported.
            Note that the ID only makes sense for AND/OR groups.
            

            C++ signature :
                unsigned int GetWriteId(RDKit::StereoGroup {lvalue})
        """
    @staticmethod
    def SetWriteId( arg1: StereoGroup, arg2: int) -> None: 
        """
        SetWriteId( arg1: StereoGroup, arg2: int) -> None
            return the StereoGroup's ID that will be exported.
            Note that the ID only makes sense for AND/OR groups.
            

            C++ signature :
                void SetWriteId(RDKit::StereoGroup {lvalue},unsigned int)
        """
    pass
class StereoGroupType(Boost.Python.enum, int):
    STEREO_ABSOLUTE = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
    STEREO_AND = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
    STEREO_OR = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
    __slots__ = ()
    names = {'STEREO_ABSOLUTE': rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE, 'STEREO_OR': rdkit.Chem.rdchem.StereoGroupType.STEREO_OR, 'STEREO_AND': rdkit.Chem.rdchem.StereoGroupType.STEREO_AND}
    values = {0: rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE, 1: rdkit.Chem.rdchem.StereoGroupType.STEREO_OR, 2: rdkit.Chem.rdchem.StereoGroupType.STEREO_AND}
    pass
class StereoGroup_vect(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: StereoGroup_vect, arg2: object) -> bool: 
        """
        __contains__( arg1: StereoGroup_vect, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: StereoGroup_vect, arg2: object) -> None: 
        """
        __delitem__( arg1: StereoGroup_vect, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<RDKit::StereoGroup*, std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > > > __iter__(boost::python::back_reference<std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> >&>)
        """
    @staticmethod
    def __len__( arg1: StereoGroup_vect) -> int: 
        """
        __len__( arg1: StereoGroup_vect) -> int

            C++ signature :
                unsigned long __len__(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: StereoGroup_vect, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: StereoGroup_vect, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: StereoGroup_vect, arg2: object) -> None: 
        """
        append( arg1: StereoGroup_vect, arg2: object) -> None

            C++ signature :
                void append(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: StereoGroup_vect, arg2: object) -> None: 
        """
        extend( arg1: StereoGroup_vect, arg2: object) -> None

            C++ signature :
                void extend(std::vector<RDKit::StereoGroup, std::allocator<RDKit::StereoGroup> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class StereoInfo(Boost.Python.instance):
    """
    Class describing stereochemistry
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def centeredOn(self) -> None:
        """
        index of the item the stereo concerns

        :type: None
        """
    @property
    def controllingAtoms(self) -> None:
        """
        indices of the atoms controlling the stereo

        :type: None
        """
    @property
    def descriptor(self) -> None:
        """
        stereo descriptor

        :type: None
        """
    @property
    def permutation(self) -> None:
        """
        permutation index (used for non-tetrahedral chirality)

        :type: None
        """
    @property
    def specified(self) -> None:
        """
        whether or not it is specified

        :type: None
        """
    @property
    def type(self) -> None:
        """
        the type of stereo

        :type: None
        """
    NOATOM = 4294967295
    __instance_size__ = 72
    pass
class StereoSpecified(Boost.Python.enum, int):
    Specified = rdkit.Chem.rdchem.StereoSpecified.Specified
    Unknown = rdkit.Chem.rdchem.StereoSpecified.Unknown
    Unspecified = rdkit.Chem.rdchem.StereoSpecified.Unspecified
    __slots__ = ()
    names = {'Unspecified': rdkit.Chem.rdchem.StereoSpecified.Unspecified, 'Specified': rdkit.Chem.rdchem.StereoSpecified.Specified, 'Unknown': rdkit.Chem.rdchem.StereoSpecified.Unknown}
    values = {0: rdkit.Chem.rdchem.StereoSpecified.Unspecified, 1: rdkit.Chem.rdchem.StereoSpecified.Specified, 2: rdkit.Chem.rdchem.StereoSpecified.Unknown}
    pass
class StereoType(Boost.Python.enum, int):
    Atom_Octahedral = rdkit.Chem.rdchem.StereoType.Atom_Octahedral
    Atom_SquarePlanar = rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar
    Atom_Tetrahedral = rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral
    Atom_TrigonalBipyramidal = rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal
    Bond_Atropisomer = rdkit.Chem.rdchem.StereoType.Bond_Atropisomer
    Bond_Cumulene_Even = rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even
    Bond_Double = rdkit.Chem.rdchem.StereoType.Bond_Double
    Unspecified = rdkit.Chem.rdchem.StereoType.Unspecified
    __slots__ = ()
    names = {'Unspecified': rdkit.Chem.rdchem.StereoType.Unspecified, 'Atom_Tetrahedral': rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral, 'Atom_SquarePlanar': rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar, 'Atom_TrigonalBipyramidal': rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal, 'Atom_Octahedral': rdkit.Chem.rdchem.StereoType.Atom_Octahedral, 'Bond_Double': rdkit.Chem.rdchem.StereoType.Bond_Double, 'Bond_Cumulene_Even': rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even, 'Bond_Atropisomer': rdkit.Chem.rdchem.StereoType.Bond_Atropisomer}
    values = {0: rdkit.Chem.rdchem.StereoType.Unspecified, 1: rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral, 2: rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar, 3: rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal, 4: rdkit.Chem.rdchem.StereoType.Atom_Octahedral, 5: rdkit.Chem.rdchem.StereoType.Bond_Double, 6: rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even, 7: rdkit.Chem.rdchem.StereoType.Bond_Atropisomer}
    pass
class SubstanceGroup(Boost.Python.instance):
    """
    A collection of atoms and bonds with associated properties
    """
    @staticmethod
    def AddAtomWithBookmark( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddAtomWithBookmark( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddAtomWithBookmark(RDKit::SubstanceGroup {lvalue},int)
        """
    @staticmethod
    def AddAtomWithIdx( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddAtomWithIdx( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddAtomWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
    @staticmethod
    def AddAttachPoint( arg1: SubstanceGroup, arg2: int, arg3: int, arg4: str) -> None: 
        """
        AddAttachPoint( arg1: SubstanceGroup, arg2: int, arg3: int, arg4: str) -> None

            C++ signature :
                void AddAttachPoint(RDKit::SubstanceGroup {lvalue},unsigned int,int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def AddBondWithBookmark( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddBondWithBookmark( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddBondWithBookmark(RDKit::SubstanceGroup {lvalue},int)
        """
    @staticmethod
    def AddBondWithIdx( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddBondWithIdx( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddBondWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
    @staticmethod
    def AddBracket( arg1: SubstanceGroup, arg2: object) -> None: 
        """
        AddBracket( arg1: SubstanceGroup, arg2: object) -> None

            C++ signature :
                void AddBracket(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
    @staticmethod
    def AddCState( arg1: SubstanceGroup, arg2: int, arg3: Point3D) -> None: 
        """
        AddCState( arg1: SubstanceGroup, arg2: int, arg3: Point3D) -> None

            C++ signature :
                void AddCState(RDKit::SubstanceGroup {lvalue},unsigned int,RDGeom::Point3D)
        """
    @staticmethod
    def AddParentAtomWithBookmark( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddParentAtomWithBookmark( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddParentAtomWithBookmark(RDKit::SubstanceGroup {lvalue},int)
        """
    @staticmethod
    def AddParentAtomWithIdx( arg1: SubstanceGroup, arg2: int) -> None: 
        """
        AddParentAtomWithIdx( arg1: SubstanceGroup, arg2: int) -> None

            C++ signature :
                void AddParentAtomWithIdx(RDKit::SubstanceGroup {lvalue},unsigned int)
        """
    @staticmethod
    def ClearAttachPoints( arg1: SubstanceGroup) -> None: 
        """
        ClearAttachPoints( arg1: SubstanceGroup) -> None

            C++ signature :
                void ClearAttachPoints(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def ClearBrackets( arg1: SubstanceGroup) -> None: 
        """
        ClearBrackets( arg1: SubstanceGroup) -> None

            C++ signature :
                void ClearBrackets(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def ClearCStates( arg1: SubstanceGroup) -> None: 
        """
        ClearCStates( arg1: SubstanceGroup) -> None

            C++ signature :
                void ClearCStates(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def ClearProp( arg1: SubstanceGroup, arg2: str) -> None: 
        """
        ClearProp( arg1: SubstanceGroup, arg2: str) -> None
            Removes a particular property (does nothing if not set).
            
            

            C++ signature :
                void ClearProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetAtoms( arg1: SubstanceGroup) -> _vectj: 
        """
        GetAtoms( arg1: SubstanceGroup) -> _vectj
            returns a list of the indices of the atoms in this SubstanceGroup

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetAtoms(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def GetAttachPoints( arg1: SubstanceGroup) -> tuple: 
        """
        GetAttachPoints( arg1: SubstanceGroup) -> tuple

            C++ signature :
                boost::python::tuple GetAttachPoints(RDKit::SubstanceGroup)
        """
    @staticmethod
    def GetBonds( arg1: SubstanceGroup) -> _vectj: 
        """
        GetBonds( arg1: SubstanceGroup) -> _vectj
            returns a list of the indices of the bonds in this SubstanceGroup

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetBonds(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def GetBoolProp( arg1: SubstanceGroup, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: SubstanceGroup, arg2: str) -> bool
            returns the value of a particular property

            C++ signature :
                bool GetBoolProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetBrackets( arg1: SubstanceGroup) -> tuple: 
        """
        GetBrackets( arg1: SubstanceGroup) -> tuple

            C++ signature :
                boost::python::tuple GetBrackets(RDKit::SubstanceGroup)
        """
    @staticmethod
    def GetCStates( arg1: SubstanceGroup) -> tuple: 
        """
        GetCStates( arg1: SubstanceGroup) -> tuple

            C++ signature :
                boost::python::tuple GetCStates(RDKit::SubstanceGroup)
        """
    @staticmethod
    def GetDoubleProp( arg1: SubstanceGroup, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: SubstanceGroup, arg2: str) -> float
            returns the value of a particular property

            C++ signature :
                double GetDoubleProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetIndexInMol( arg1: SubstanceGroup) -> int: 
        """
        GetIndexInMol( arg1: SubstanceGroup) -> int
            returns the index of this SubstanceGroup in the owning molecule's list.

            C++ signature :
                unsigned int GetIndexInMol(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def GetIntProp( arg1: SubstanceGroup, arg2: str) -> int: 
        """
        GetIntProp( arg1: SubstanceGroup, arg2: str) -> int
            returns the value of a particular property

            C++ signature :
                int GetIntProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetOwningMol( arg1: SubstanceGroup) -> Mol: 
        """
        GetOwningMol( arg1: SubstanceGroup) -> Mol
            returns the molecule owning this SubstanceGroup

            C++ signature :
                RDKit::ROMol {lvalue} GetOwningMol(RDKit::SubstanceGroup {lvalue})
        """
    @staticmethod
    def GetParentAtoms( arg1: SubstanceGroup) -> _vectj: 
        """
        GetParentAtoms( arg1: SubstanceGroup) -> _vectj
            returns a list of the indices of the parent atoms in this SubstanceGroup

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetParentAtoms(RDKit::SubstanceGroup {lvalue})
        """
    def GetProp(self, key: str, autoConvert: bool = False) -> object: 
        """
        GetProp( self: SubstanceGroup, key: str, autoConvert: bool = False) -> object
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                boost::python::api::object GetProp(RDKit::SubstanceGroup const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: SubstanceGroup, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a list of the properties set on the SubstanceGroup.
            
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::SubstanceGroup {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: SubstanceGroup, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict
            Returns a dictionary of the properties set on the SubstanceGroup.
             n.b. some properties cannot be converted to python types.
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::SubstanceGroup [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    def GetStringVectProp( arg1: SubstanceGroup, arg2: str) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetStringVectProp( arg1: SubstanceGroup, arg2: str) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            returns the value of a particular property

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetStringVectProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetUnsignedProp( arg1: SubstanceGroup, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: SubstanceGroup, arg2: str) -> int
            returns the value of a particular property

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def GetUnsignedVectProp( arg1: SubstanceGroup, arg2: str) -> _vectj: 
        """
        GetUnsignedVectProp( arg1: SubstanceGroup, arg2: str) -> _vectj
            returns the value of a particular property

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetUnsignedVectProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def HasProp( arg1: SubstanceGroup, arg2: str) -> bool: 
        """
        HasProp( arg1: SubstanceGroup, arg2: str) -> bool
            returns whether or not a particular property exists

            C++ signature :
                bool HasProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def SetAtoms( arg1: SubstanceGroup, arg2: object) -> None: 
        """
        SetAtoms( arg1: SubstanceGroup, arg2: object) -> None
            Set the list of the indices of the atoms in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.

            C++ signature :
                void SetAtoms(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
    @staticmethod
    def SetBonds( arg1: SubstanceGroup, arg2: object) -> None: 
        """
        SetBonds( arg1: SubstanceGroup, arg2: object) -> None
            Set the list of the indices of the bonds in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.

            C++ signature :
                void SetBonds(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None: 
        """
        SetBoolProp( self: SubstanceGroup, key: str, val: bool, computed: bool = False) -> None
            sets the value of a particular property

            C++ signature :
                void SetBoolProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None: 
        """
        SetDoubleProp( self: SubstanceGroup, key: str, val: float, computed: bool = False) -> None
            sets the value of a particular property

            C++ signature :
                void SetDoubleProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False])
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetIntProp( self: SubstanceGroup, key: str, val: int, computed: bool = False) -> None
            sets the value of a particular property

            C++ signature :
                void SetIntProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False])
        """
    @staticmethod
    def SetParentAtoms( arg1: SubstanceGroup, arg2: object) -> None: 
        """
        SetParentAtoms( arg1: SubstanceGroup, arg2: object) -> None
            Set the list of the indices of the parent atoms in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.

            C++ signature :
                void SetParentAtoms(RDKit::SubstanceGroup {lvalue},boost::python::api::object)
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None: 
        """
        SetProp( self: SubstanceGroup, key: str, val: str, computed: bool = False) -> None
            sets the value of a particular property

            C++ signature :
                void SetProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetUnsignedProp( self: SubstanceGroup, key: str, val: int, computed: bool = False) -> None
            sets the value of a particular property

            C++ signature :
                void SetUnsignedProp(RDKit::SubstanceGroup {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int [,bool=False])
        """
    pass
class SubstanceGroupAttach(Boost.Python.instance):
    """
    AttachPoint for a SubstanceGroup
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def aIdx(self) -> None:
        """
        attachment index

        :type: None
        """
    @property
    def id(self) -> None:
        """
        attachment id

        :type: None
        """
    @property
    def lvIdx(self) -> None:
        """
        leaving atom or index (0 for implied)

        :type: None
        """
    __instance_size__ = 40
    pass
class SubstanceGroupCState(Boost.Python.instance):
    """
    CSTATE for a SubstanceGroup
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def bondIdx(self) -> None:
        """
        :type: None
        """
    @property
    def vector(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 40
    pass
class SubstanceGroup_VECT(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: SubstanceGroup_VECT, arg2: object) -> bool: 
        """
        __contains__( arg1: SubstanceGroup_VECT, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: SubstanceGroup_VECT, arg2: object) -> None: 
        """
        __delitem__( arg1: SubstanceGroup_VECT, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<RDKit::SubstanceGroup*, std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > > > __iter__(boost::python::back_reference<std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> >&>)
        """
    @staticmethod
    def __len__( arg1: SubstanceGroup_VECT) -> int: 
        """
        __len__( arg1: SubstanceGroup_VECT) -> int

            C++ signature :
                unsigned long __len__(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: SubstanceGroup_VECT, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: SubstanceGroup_VECT, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: SubstanceGroup_VECT, arg2: object) -> None: 
        """
        append( arg1: SubstanceGroup_VECT, arg2: object) -> None

            C++ signature :
                void append(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: SubstanceGroup_VECT, arg2: object) -> None: 
        """
        extend( arg1: SubstanceGroup_VECT, arg2: object) -> None

            C++ signature :
                void extend(std::vector<RDKit::SubstanceGroup, std::allocator<RDKit::SubstanceGroup> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class SubstructMatchParameters(Boost.Python.instance):
    """
    Parameters controlling substructure matching
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def setExtraFinalCheck( arg1: SubstructMatchParameters, arg2: object) -> None: 
        """
        setExtraFinalCheck( arg1: SubstructMatchParameters, arg2: object) -> None
            allows you to provide a function that will be called
                           with the molecule
                       and a vector of atom IDs containing a potential match.
                       The function should return true or false indicating whether or not
                       that match should be accepted.

            C++ signature :
                void setExtraFinalCheck(RDKit::SubstructMatchParameters {lvalue},boost::python::api::object)
        """
    @property
    def aromaticMatchesConjugated(self) -> None:
        """
        aromatic and conjugated bonds match each other

        :type: None
        """
    @property
    def atomProperties(self) -> None:
        """
        atom properties that must be equivalent in order to match.

        :type: None
        """
    @property
    def bondProperties(self) -> None:
        """
        bond properties that must be equivalent in order to match.

        :type: None
        """
    @property
    def maxMatches(self) -> None:
        """
        maximum number of matches to return

        :type: None
        """
    @property
    def numThreads(self) -> None:
        """
        number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware.

        :type: None
        """
    @property
    def recursionPossible(self) -> None:
        """
        Allow recursive queries

        :type: None
        """
    @property
    def uniquify(self) -> None:
        """
        :type: None
        """
    @property
    def useChirality(self) -> None:
        """
        Use chirality in determining whether or not atoms/bonds match

        :type: None
        """
    @property
    def useEnhancedStereo(self) -> None:
        """
        take enhanced stereochemistry into account while doing the match. This only has an effect if useChirality is also True.

        :type: None
        """
    @property
    def useGenericMatchers(self) -> None:
        """
        use generic groups (=homology groups) as a post-filtering step (if any are present in the molecule)

        :type: None
        """
    @property
    def useQueryQueryMatches(self) -> None:
        """
        Consider query-query matches, not just simple matches

        :type: None
        """
    __instance_size__ = 120
    pass
class _ROConformerSeq(Boost.Python.instance):
    """
    Read-only sequence of conformers, not constructible from Python.
    """
    @staticmethod
    def __getitem__( arg1: _ROConformerSeq, arg2: int) -> Conformer: 
        """
        __getitem__( arg1: _ROConformerSeq, arg2: int) -> Conformer

            C++ signature :
                RDKit::Conformer* __getitem__(RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor> {lvalue},int)
        """
    @staticmethod
    def __iter__( arg1: _ROConformerSeq) -> _ROConformerSeq: 
        """
        __iter__( arg1: _ROConformerSeq) -> _ROConformerSeq

            C++ signature :
                RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor>* __iter__(RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor> {lvalue})
        """
    @staticmethod
    def __len__( arg1: _ROConformerSeq) -> int: 
        """
        __len__( arg1: _ROConformerSeq) -> int

            C++ signature :
                int __len__(RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor> {lvalue})
        """
    @staticmethod
    def __next__( arg1: _ROConformerSeq) -> Conformer: 
        """
        __next__( arg1: _ROConformerSeq) -> Conformer

            C++ signature :
                RDKit::Conformer* __next__(RDKit::ReadOnlySeq<std::_List_iterator<boost::shared_ptr<RDKit::Conformer> >, boost::shared_ptr<RDKit::Conformer>&, RDKit::ConformerCountFunctor> {lvalue})
        """
    pass
class _ROQAtomSeq(Boost.Python.instance):
    """
    Read-only sequence of atoms matching a query, not constructible from Python.
    """
    @staticmethod
    def __getitem__( arg1: _ROQAtomSeq, arg2: int) -> Atom: 
        """
        __getitem__( arg1: _ROQAtomSeq, arg2: int) -> Atom

            C++ signature :
                RDKit::Atom* __getitem__(RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor> {lvalue},int)
        """
    @staticmethod
    def __iter__( arg1: _ROQAtomSeq) -> _ROQAtomSeq: 
        """
        __iter__( arg1: _ROQAtomSeq) -> _ROQAtomSeq

            C++ signature :
                RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor>* __iter__(RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor> {lvalue})
        """
    @staticmethod
    def __len__( arg1: _ROQAtomSeq) -> int: 
        """
        __len__( arg1: _ROQAtomSeq) -> int

            C++ signature :
                int __len__(RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor> {lvalue})
        """
    @staticmethod
    def __next__( arg1: _ROQAtomSeq) -> Atom: 
        """
        __next__( arg1: _ROQAtomSeq) -> Atom

            C++ signature :
                RDKit::Atom* __next__(RDKit::ReadOnlySeq<RDKit::QueryAtomIterator_<RDKit::Atom, RDKit::ROMol>, RDKit::Atom*, RDKit::AtomCountFunctor> {lvalue})
        """
    pass
class _cppMolSanitizeException(Boost.Python.instance):
    """
    exception arising from sanitization
    """
    @staticmethod
    def GetType( arg1: _cppMolSanitizeException) -> str: 
        """
        GetType( arg1: _cppMolSanitizeException) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetType(RDKit::MolSanitizeException {lvalue})
        """
    @staticmethod
    def Message( arg1: _cppMolSanitizeException) -> str: 
        """
        Message( arg1: _cppMolSanitizeException) -> str

            C++ signature :
                char const* Message(RDKit::MolSanitizeException {lvalue})
        """
    pass
class _cppAtomSanitizeException(_cppMolSanitizeException, Boost.Python.instance):
    """
    exception arising from sanitization
    """
    @staticmethod
    def GetAtomIdx( arg1: _cppAtomSanitizeException) -> int: 
        """
        GetAtomIdx( arg1: _cppAtomSanitizeException) -> int

            C++ signature :
                unsigned int GetAtomIdx(RDKit::AtomSanitizeException {lvalue})
        """
    pass
class _cppAtomValenceException(_cppAtomSanitizeException, _cppMolSanitizeException, Boost.Python.instance):
    """
    exception arising from sanitization
    """
    pass
class _cppAtomKekulizeException(_cppMolSanitizeException, Boost.Python.instance):
    """
    exception arising from sanitization
    """
    @staticmethod
    def GetAtomIndices( arg1: _cppAtomKekulizeException) -> tuple: 
        """
        GetAtomIndices( arg1: _cppAtomKekulizeException) -> tuple

            C++ signature :
                boost::python::tuple GetAtomIndices(RDKit::KekulizeException)
        """
    pass
class _listN5boost10shared_ptrIN5RDKit9ConformerEEE(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object) -> bool: 
        """
        __contains__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object) -> None: 
        """
        __delitem__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object) -> None

            C++ signature :
                void __delitem__(std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::_List_iterator<boost::shared_ptr<RDKit::Conformer> > > __iter__(boost::python::back_reference<std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > >&>)
        """
    @staticmethod
    def __len__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE) -> int: 
        """
        __len__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE) -> int

            C++ signature :
                unsigned long __len__(std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: _listN5boost10shared_ptrIN5RDKit9ConformerEEE, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::__cxx11::list<boost::shared_ptr<RDKit::Conformer>, std::allocator<boost::shared_ptr<RDKit::Conformer> > > {lvalue},_object*,_object*)
        """
    __instance_size__ = 48
    pass
class _listPN5RDKit4AtomE(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: _listPN5RDKit4AtomE, arg2: object) -> bool: 
        """
        __contains__( arg1: _listPN5RDKit4AtomE, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: _listPN5RDKit4AtomE, arg2: object) -> None: 
        """
        __delitem__( arg1: _listPN5RDKit4AtomE, arg2: object) -> None

            C++ signature :
                void __delitem__(std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::_List_iterator<RDKit::Atom*> > __iter__(boost::python::back_reference<std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> >&>)
        """
    @staticmethod
    def __len__( arg1: _listPN5RDKit4AtomE) -> int: 
        """
        __len__( arg1: _listPN5RDKit4AtomE) -> int

            C++ signature :
                unsigned long __len__(std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: _listPN5RDKit4AtomE, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: _listPN5RDKit4AtomE, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::__cxx11::list<RDKit::Atom*, std::allocator<RDKit::Atom*> > {lvalue},_object*,_object*)
        """
    __instance_size__ = 48
    pass
class _listPN5RDKit4BondE(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: _listPN5RDKit4BondE, arg2: object) -> bool: 
        """
        __contains__( arg1: _listPN5RDKit4BondE, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: _listPN5RDKit4BondE, arg2: object) -> None: 
        """
        __delitem__( arg1: _listPN5RDKit4BondE, arg2: object) -> None

            C++ signature :
                void __delitem__(std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::_List_iterator<RDKit::Bond*> > __iter__(boost::python::back_reference<std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> >&>)
        """
    @staticmethod
    def __len__( arg1: _listPN5RDKit4BondE) -> int: 
        """
        __len__( arg1: _listPN5RDKit4BondE) -> int

            C++ signature :
                unsigned long __len__(std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: _listPN5RDKit4BondE, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: _listPN5RDKit4BondE, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::__cxx11::list<RDKit::Bond*, std::allocator<RDKit::Bond*> > {lvalue},_object*,_object*)
        """
    __instance_size__ = 48
    pass
def AddMolSubstanceGroup( mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup:
    """
    AddMolSubstanceGroup( mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup
        adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup

        C++ signature :
            RDKit::SubstanceGroup* AddMolSubstanceGroup(RDKit::ROMol {lvalue},RDKit::SubstanceGroup)
    """
def ClearMolSubstanceGroups( arg1: Mol) -> None:
    """
    ClearMolSubstanceGroups( arg1: Mol) -> None
        removes all SubstanceGroups from a molecule (if any)

        C++ signature :
            void ClearMolSubstanceGroups(RDKit::ROMol {lvalue})
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
def ForwardStereoGroupIds( arg1: Mol) -> None:
    """
    ForwardStereoGroupIds( arg1: Mol) -> None
        Forward the original Stereo Group IDs when exporting the Mol.

        C++ signature :
            void ForwardStereoGroupIds(RDKit::ROMol {lvalue})
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
def GetPeriodicTable() -> PeriodicTable:
    """
    GetPeriodicTable() -> PeriodicTable
        Returns the application's PeriodicTable instance.
        
        

        C++ signature :
            RDKit::PeriodicTable* GetPeriodicTable()
    """
def GetSupplementalSmilesLabel( atom: Atom) -> str:
    """
    GetSupplementalSmilesLabel( atom: Atom) -> str
        Gets the supplemental smiles label on an atom, returns an empty string if not present.

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetSupplementalSmilesLabel(RDKit::Atom const*)
    """
def MolBundleCanSerialize() -> bool:
    """
    MolBundleCanSerialize() -> bool
        Returns True if the MolBundle is serializable (requires boost serialization

        C++ signature :
            bool MolBundleCanSerialize()
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
def SetDefaultPickleProperties( arg1: int) -> None:
    """
    SetDefaultPickleProperties( arg1: int) -> None
        Set the current global mol pickler options.

        C++ signature :
            void SetDefaultPickleProperties(unsigned int)
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
def tossit() -> None:
    """
    tossit() -> None

        C++ signature :
            void tossit()
    """
ALLOW_CHARGE_SEPARATION = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
ALLOW_INCOMPLETE_OCTETS = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
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
KEKULE_ALL = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
MolProps = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
NoConformers = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
STEREO_ABSOLUTE = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
UNCONSTRAINED_ANIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
