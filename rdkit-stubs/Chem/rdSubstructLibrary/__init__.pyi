from __future__ import annotations
import rdkit.Chem.rdSubstructLibrary
import typing
import Boost.Python

__all__ = [
    "AddPatterns",
    "CachedMolHolder",
    "CachedSmilesMolHolder",
    "CachedTrustedSmilesMolHolder",
    "FPHolderBase",
    "KeyFromPropHolder",
    "KeyHolderBase",
    "MolHolder",
    "MolHolderBase",
    "PatternHolder",
    "SubstructLibrary",
    "SubstructLibraryCanSerialize",
    "TautomerPatternHolder"
]


class MolHolderBase(Boost.Python.instance):
    @staticmethod
    def AddMol( arg1: MolHolderBase, arg2: Mol) -> int: 
        """
        AddMol( arg1: MolHolderBase, arg2: Mol) -> int
            Adds molecule to the molecule holder

            C++ signature :
                unsigned int AddMol(RDKit::MolHolderBase {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def GetMol( arg1: MolHolderBase, arg2: int) -> Mol: 
        """
        GetMol( arg1: MolHolderBase, arg2: int) -> Mol
            Returns a particular molecule in the molecule holder
            
              ARGUMENTS:
                - idx: which molecule to return
            
                - sanitize: if sanitize is False, return the internal molecule state [default True]
            
              NOTE: molecule indices start at 0
            

            C++ signature :
                boost::shared_ptr<RDKit::ROMol> GetMol(RDKit::MolHolderBase {lvalue},unsigned int)
        """
    @staticmethod
    def __len__( arg1: MolHolderBase) -> int: 
        """
        __len__( arg1: MolHolderBase) -> int

            C++ signature :
                unsigned int __len__(RDKit::MolHolderBase {lvalue})

            C++ signature :
                unsigned int __len__(RDKit::MolHolderBase {lvalue})
        """
    pass
class CachedSmilesMolHolder(MolHolderBase, Boost.Python.instance):
    """
    Holds molecules as smiles string
    This allows more molecules to be held in memory at a time
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                           The smiles is stored as-is, no checking is done for validity.
      GetMol(idx) -> return the molecule at index idx
    """
    @staticmethod
    def AddSmiles( arg1: CachedSmilesMolHolder, smiles: str) -> int: 
        """
        AddSmiles( arg1: CachedSmilesMolHolder, smiles: str) -> int
            Add a trusted smiles string to the molecule holder, no checking is done on the input data

            C++ signature :
                unsigned int AddSmiles(RDKit::CachedSmilesMolHolder {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 40
    pass
class CachedTrustedSmilesMolHolder(MolHolderBase, Boost.Python.instance):
    """
    Holds molecules as trusted smiles string
    This allows more molecules to be held in memory at a time and avoids RDKit sanitization
    overhead.
    See: http://rdkit.blogspot.com/2016/09/avoiding-unnecessary-work-and.html
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddSmiles(smiles) -> adds a smiles string to the molecule holder, returns index of molecule
                           The smiles is stored as-is, no checking is done for validity.
      GetMol(idx,s) -> return the molecule at index idx, 
                  note, only light sanitization is done here, for instance
                  the molecules RingInfo is not initialized
    """
    @staticmethod
    def AddSmiles( arg1: CachedTrustedSmilesMolHolder, smiles: str) -> int: 
        """
        AddSmiles( arg1: CachedTrustedSmilesMolHolder, smiles: str) -> int
            Add a trusted smiles string to the molecule holder, no checking is done on the input data

            C++ signature :
                unsigned int AddSmiles(RDKit::CachedTrustedSmilesMolHolder {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 40
    pass
class FPHolderBase(Boost.Python.instance):
    @staticmethod
    def AddFingerprint( arg1: FPHolderBase, arg2: ExplicitBitVect) -> int: 
        """
        AddFingerprint( arg1: FPHolderBase, arg2: ExplicitBitVect) -> int
            Adds a raw bit vector to the fingerprint database, returns the index of the supplied pattern

            C++ signature :
                unsigned int AddFingerprint(RDKit::FPHolderBase {lvalue},ExplicitBitVect)
        """
    @staticmethod
    def AddMol( arg1: FPHolderBase, arg2: Mol) -> int: 
        """
        AddMol( arg1: FPHolderBase, arg2: Mol) -> int
            Adds a molecule to the fingerprint database, returns the index of the new pattern

            C++ signature :
                unsigned int AddMol(RDKit::FPHolderBase {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def GetFingerprint( arg1: FPHolderBase, arg2: int) -> ExplicitBitVect: 
        """
        GetFingerprint( arg1: FPHolderBase, arg2: int) -> ExplicitBitVect
            Return the bit vector at the specified index

            C++ signature :
                ExplicitBitVect GetFingerprint(RDKit::FPHolderBase {lvalue},unsigned int)
        """
    @staticmethod
    def MakeFingerprint( arg1: FPHolderBase, mol: Mol) -> ExplicitBitVect: 
        """
        MakeFingerprint( arg1: FPHolderBase, mol: Mol) -> ExplicitBitVect
            Compute the query bits for the holder

            C++ signature :
                ExplicitBitVect* MakeFingerprint(RDKit::FPHolderBase {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def PassesFilter( arg1: FPHolderBase, idx: int, query: ExplicitBitVect) -> bool: 
        """
        PassesFilter( arg1: FPHolderBase, idx: int, query: ExplicitBitVect) -> bool
            Returns True if the specified index passes the filter supplied by the query bit vector

            C++ signature :
                bool PassesFilter(RDKit::FPHolderBase {lvalue},unsigned int,ExplicitBitVect)
        """
    @staticmethod
    def __len__( arg1: FPHolderBase) -> int: 
        """
        __len__( arg1: FPHolderBase) -> int

            C++ signature :
                unsigned int __len__(RDKit::FPHolderBase {lvalue})
        """
    pass
class KeyHolderBase(Boost.Python.instance):
    @staticmethod
    def AddKey( arg1: KeyHolderBase, arg2: str) -> int: 
        """
        AddKey( arg1: KeyHolderBase, arg2: str) -> int
            Add a key to the key holder, must be manually synced

            C++ signature :
                unsigned int AddKey(RDKit::KeyHolderBase {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def AddMol( arg1: KeyHolderBase, arg2: Mol) -> int: 
        """
        AddMol( arg1: KeyHolderBase, arg2: Mol) -> int
            Adds a molecule to the fingerprint database, returns the index of the new pattern

            C++ signature :
                unsigned int AddMol(RDKit::KeyHolderBase {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def GetKey( arg1: KeyHolderBase, arg2: int) -> str: 
        """
        GetKey( arg1: KeyHolderBase, arg2: int) -> str
            Return the key at the specified index

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetKey(RDKit::KeyHolderBase {lvalue},unsigned int)
        """
    @staticmethod
    def GetKeys( arg1: KeyHolderBase, arg2: _vectj) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetKeys( arg1: KeyHolderBase, arg2: _vectj) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns the keys for the given indices as return by GetMatches 
            
              ARGUMENTS:
                - indices: The indices of the keys
            
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetKeys(RDKit::KeyHolderBase {lvalue},std::vector<unsigned int, std::allocator<unsigned int> >)
        """
    @staticmethod
    def __len__( arg1: KeyHolderBase) -> int: 
        """
        __len__( arg1: KeyHolderBase) -> int

            C++ signature :
                unsigned int __len__(RDKit::KeyHolderBase {lvalue})
        """
    pass
class KeyFromPropHolder(KeyHolderBase, Boost.Python.instance):
    """
    Holds keys to return external references to the molecules in the molholder.
    By default use the _Name property but can be overridden to be any property
    """
    @staticmethod
    def GetPropName( arg1: KeyFromPropHolder) -> str: 
        """
        GetPropName( arg1: KeyFromPropHolder) -> str
            Return the key for the given molecule index

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetPropName(RDKit::KeyFromPropHolder {lvalue})
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
    def __init__( arg1: object, arg2: str) -> None: ...
    __instance_size__ = 40
    pass
class MolHolder(MolHolderBase, Boost.Python.instance):
    """
    Holds raw in-memory molecules
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule
      GetMol(idx,sanitize=True) -> return the molecule at index idx
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 40
    pass
class CachedMolHolder(MolHolderBase, Boost.Python.instance):
    """
    Holds molecules in their binary representation.
    This allows more molecules to be held in memory at a time
      AddMol(mol) -> adds a molecule to the molecule holder, returns index of molecule

      AddBinary(data) -> adds a picked molecule molecule to the molecule holder, returns index of molecule
                         The data is stored as-is, no checking is done for validity.
      GetMol(idx) -> return the molecule at index idx
    """
    @staticmethod
    def AddBinary( arg1: CachedMolHolder, pickle: str) -> int: 
        """
        AddBinary( arg1: CachedMolHolder, pickle: str) -> int
            Add a binary pickle to the molecule holder, no checking is done on the input data

            C++ signature :
                unsigned int AddBinary(RDKit::CachedMolHolder {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 40
    pass
class PatternHolder(FPHolderBase, Boost.Python.instance):
    """
    Holds fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
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
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: ...
    __instance_size__ = 40
    pass
class SubstructLibrary(Boost.Python.instance):
    """
    SubstructLibrary: This provides a simple API for substructure searching large datasets
    The SubstructLibrary takes full advantage of available threads during the search operation.
    Basic operation is simple

    >>> from __future__ import print_function
    >>> import os
    >>> from rdkit import Chem, RDConfig
    >>> from rdkit.Chem import rdSubstructLibrary
    >>> library = rdSubstructLibrary.SubstructLibrary()
    >>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir, 
    ...                               'NCI', 'first_200.props.sdf')):
    ...   idx = library.AddMol(mol)
    >>> core = Chem.MolFromSmarts('CCCCOC')
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    Substructure matching options can be sent into GetMatches:

    >>> indices = library.GetMatches(core, useChirality=False) 
    >>> len(indices)
    11

    Controlling the number of threads or the maximum number of matches returned:
    is also available (the default is to run on all cores)

    >>> indices = library.GetMatches(core, numThreads=2, maxResults=10) 
    >>> len(indices)
    10

    Working on larger datasets:

    Molecules are fairly large objects and will limit the number that can be kept in memory.
    To assist this we supply three other molecule holders:
      CachedMolHolder - stores molecules as their pickled representation

      CachedSmilesMolHolder - stores molecules internally as smiles strings

      CachedTrustedSmilesMolHolder = excepts (and stores) molecules as trusted smiles strings

    Using Pattern fingerprints as a pre-filter:
    Pattern fingerprints provide an easy way to indicate whether the substructure search should be
    be done at all.  This is particularly useful with the Binary and Smiles based molecule holders
    as they have an expensive molecule creation step in addition to the substructure searching step
     
    >>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.CachedSmilesMolHolder(), 
    ...                                               rdSubstructLibrary.PatternHolder())
    >>> for mol in Chem.SDMolSupplier(os.path.join(RDConfig.RDDataDir, 
    ...                               'NCI', 'first_200.props.sdf')):
    ...   idx = library.AddMol(mol)
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    This (obviously) takes longer to initialize.  However, both the molecule and pattern
    holders can be populated with raw data, a simple example is below:

    >>> import csv
    >>> molholder = rdSubstructLibrary.CachedSmilesMolHolder()
    >>> pattern_holder = rdSubstructLibrary.PatternHolder()
    >>> with open(os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.tpsa.csv')) as inf:
    ...   for i, row in enumerate(csv.reader(inf)):
    ...     if i:
    ...       idx = molholder.AddSmiles(row[0])
    ...       idx2 = pattern_holder.AddFingerprint(
    ...           pattern_holder.MakeFingerprint(Chem.MolFromSmiles(row[0])))
    ...       assert idx==idx2
    >>> library = rdSubstructLibrary.SubstructLibrary(molholder,pattern_holder)
    >>> indices = library.GetMatches(core)
    >>> len(indices)
    11

    Finally, the KeyFromPropHolder can be used to use external keys such as
    compound names.  By default the holder uses the '_Name' property but can
    be changed to any property.

    >>> library = rdSubstructLibrary.SubstructLibrary(rdSubstructLibrary.MolHolder(), rdSubstructLibrary.KeyFromPropHolder())
    >>> m = Chem.MolFromSmiles('CCC')
    >>> m.SetProp('_Name', 'Z11234')
    >>> idx = library.AddMol(m)
    >>> indices = library.GetMatches(m)
    >>> list(library.GetKeyHolder().GetKeys(indices))
    ['Z11234']
    """
    @staticmethod
    def AddMol( arg1: SubstructLibrary, mol: Mol) -> int: 
        """
        AddMol( arg1: SubstructLibrary, mol: Mol) -> int
            Adds a molecule to the substruct library

            C++ signature :
                unsigned int AddMol(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol)
        """
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: 
        """
        CountMatches( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int
            Get the matches for the query.
            
             Arguments:
              - query:      substructure query
              - numThreads: number of threads to use, -1 means all threads
            

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                unsigned int CountMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])
        """
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: Mol, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: object, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: object, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: MolBundle, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    @typing.overload
    def CountMatches( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> int: ...
    @staticmethod
    def GetFpHolder( arg1: SubstructLibrary) -> FPHolderBase: 
        """
        GetFpHolder( arg1: SubstructLibrary) -> FPHolderBase

            C++ signature :
                boost::shared_ptr<RDKit::FPHolderBase> GetFpHolder(RDKit::SubstructLibraryWrap {lvalue})
        """
    @staticmethod
    def GetKeyHolder( arg1: SubstructLibrary) -> KeyHolderBase: 
        """
        GetKeyHolder( arg1: SubstructLibrary) -> KeyHolderBase

            C++ signature :
                boost::shared_ptr<RDKit::KeyHolderBase> GetKeyHolder(RDKit::SubstructLibraryWrap {lvalue})
        """
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: 
        """
        GetMatches( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj
            Get the matches for the query.
            
             Arguments:
              - query:      substructure query
              - numThreads: number of threads to use, -1 means all threads
              - maxResults: maximum number of results to return

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1 [,int=1000]]]]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > GetMatches(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1 [,int=1000]])
        """
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: Mol, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: object, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: object, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: MolBundle, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    @typing.overload
    def GetMatches( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1, maxResults: int = 1000) -> _vectj: ...
    @staticmethod
    def GetMol( arg1: SubstructLibrary, arg2: int) -> Mol: 
        """
        GetMol( arg1: SubstructLibrary, arg2: int) -> Mol
            Returns a particular molecule in the molecule holder
            
              ARGUMENTS:
                - idx: which molecule to return
            
              NOTE: molecule indices start at 0
            

            C++ signature :
                boost::shared_ptr<RDKit::ROMol> GetMol(RDKit::SubstructLibraryWrap {lvalue},unsigned int)
        """
    @staticmethod
    def GetMolHolder( arg1: SubstructLibrary) -> MolHolderBase: 
        """
        GetMolHolder( arg1: SubstructLibrary) -> MolHolderBase

            C++ signature :
                boost::shared_ptr<RDKit::MolHolderBase> GetMolHolder(RDKit::SubstructLibraryWrap {lvalue})
        """
    @staticmethod
    def GetSearchOrder( arg1: SubstructLibrary) -> tuple: 
        """
        GetSearchOrder( arg1: SubstructLibrary) -> tuple
            Returns the search order for the library
            
              NOTE: molecule indices start at 0
            

            C++ signature :
                boost::python::tuple GetSearchOrder(RDKit::SubstructLibraryWrap)
        """
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: 
        """
        HasMatch( arg1: SubstructLibrary, query: Mol, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool
            Get the matches for the query.
            
             Arguments:
              - query:      substructure query
              - numThreads: number of threads to use, -1 means all threads
            

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::ROMol,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::TautomerQuery,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int [,bool=True [,bool=True [,bool=False [,int=-1]]]])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,RDKit::SubstructMatchParameters [,int=-1])

            C++ signature :
                bool HasMatch(RDKit::SubstructLibraryWrap {lvalue},RDKit::MolBundle,unsigned int,unsigned int,RDKit::SubstructMatchParameters [,int=-1])
        """
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: Mol, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: Mol, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: object, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: object, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: object, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: MolBundle, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, recursionPossible: bool = True, useChirality: bool = True, useQueryQueryMatches: bool = False, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: MolBundle, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    @typing.overload
    def HasMatch( arg1: SubstructLibrary, query: MolBundle, startIdx: int, endIdx: int, parameters: SubstructMatchParameters, numThreads: int = -1) -> bool: ...
    @staticmethod
    def InitFromStream( arg1: SubstructLibrary, stream: AtomPairsParameters) -> None: 
        """
        InitFromStream( arg1: SubstructLibrary, stream: AtomPairsParameters) -> None
            Deserialize a substructure library from a python bytes stream.
            Python doesn't allow seeking operations inside a unicode or string stream anymore
            so this requires opening a file in binary mode or using an io.ByteIO type object
            
              ARGUMENTS:
                - stream: a binary stream like object
            
              SubstructLibrary.Serialize already writes a binary stream
            
              >>> from rdkit.Chem import rdSubstructLibrary
              >>> import io
              >>> lib = rdSubstructLibrary.SubstructLibrary()
              >>> stream = io.BytesIO( lib.Serialize() )
              >>> lib.InitFromStream(stream)
            
               remember to write to text and read from a binary stream
              >>> with open('rdkit.sslib', 'w') as f: lib.ToStream(f)
              >>> with open('rdkit.sslib', 'rb') as f: lib.InitFromStream(f)
            

            C++ signature :
                void InitFromStream(RDKit::SubstructLibraryWrap {lvalue},boost::python::api::object {lvalue})
        """
    @staticmethod
    def Serialize( arg1: SubstructLibrary) -> object: 
        """
        Serialize( arg1: SubstructLibrary) -> object

            C++ signature :
                boost::python::api::object Serialize(RDKit::SubstructLibraryWrap)
        """
    @staticmethod
    def SetSearchOrder( arg1: SubstructLibrary, arg2: AtomPairsParameters) -> None: 
        """
        SetSearchOrder( arg1: SubstructLibrary, arg2: AtomPairsParameters) -> None
            Sets the search order for the library
            
              ARGUMENTS:
                - order: sequence of molecule indices
            
              NOTE: molecule indices start at 0
            

            C++ signature :
                void SetSearchOrder(RDKit::SubstructLibraryWrap {lvalue},boost::python::api::object)
        """
    @staticmethod
    def ToStream( arg1: SubstructLibrary, stream: AtomPairsParameters) -> None: 
        """
        ToStream( arg1: SubstructLibrary, stream: AtomPairsParameters) -> None
            Serialize a substructure library to a python text stream.
            The stream can be a file in text mode or an io.StringIO type object
            
              ARGUMENTS:
                - stream: a text or text stream like object
            
              >>> from rdkit.Chem import rdSubstructLibrary
              >>> import io
              >>> lib = rdSubstructLibrary.SubstructLibrary()
              >>> stream = io.StringIO()
              >>> lib.ToStream(stream)
            
               or
              >>> with open('rdkit.sslib', 'w') as stream:
              ...  lib.ToStream(stream)
            

            C++ signature :
                void ToStream(RDKit::SubstructLibraryWrap,boost::python::api::object {lvalue})
        """
    @staticmethod
    def __getinitargs__( arg1: SubstructLibrary) -> tuple: 
        """
        __getinitargs__( arg1: SubstructLibrary) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::SubstructLibraryWrap)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

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
                void __init__(_object*,boost::shared_ptr<RDKit::MolHolderBase>)

            C++ signature :
                void __init__(_object*,boost::shared_ptr<RDKit::MolHolderBase>,boost::shared_ptr<RDKit::FPHolderBase>)

            C++ signature :
                void __init__(_object*,boost::shared_ptr<RDKit::MolHolderBase>,boost::shared_ptr<RDKit::KeyHolderBase>)

            C++ signature :
                void __init__(_object*,boost::shared_ptr<RDKit::MolHolderBase>,boost::shared_ptr<RDKit::FPHolderBase>,boost::shared_ptr<RDKit::KeyHolderBase>)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: MolHolderBase) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: MolHolderBase, arg3: FPHolderBase) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: MolHolderBase, arg3: KeyHolderBase) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: MolHolderBase, arg3: FPHolderBase, arg4: KeyHolderBase) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    def __len__( arg1: SubstructLibrary) -> int: 
        """
        __len__( arg1: SubstructLibrary) -> int

            C++ signature :
                unsigned int __len__(RDKit::SubstructLibraryWrap {lvalue})
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 32
    __safe_for_unpickling__ = True
    pass
class TautomerPatternHolder(FPHolderBase, Boost.Python.instance):
    """
    Holds tautomeric fingerprints with optional, user-defined number of bits (default: 2048) used for filtering of molecules.
    These fingerprints are designed to be used with TautomerQueries.
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
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: int) -> None: ...
    __instance_size__ = 40
    pass
@typing.overload
def AddPatterns( sslib: SubstructLibrary, numThreads: int = 1) -> None:
    """
    AddPatterns( sslib: SubstructLibrary, numThreads: int = 1) -> None
        Add pattern fingerprints to the given library, use numThreads=-1 to use all available cores

        C++ signature :
            void AddPatterns(RDKit::SubstructLibraryWrap {lvalue} [,int=1])

        C++ signature :
            void AddPatterns(RDKit::SubstructLibraryWrap {lvalue},boost::shared_ptr<RDKit::FPHolderBase> [,int=1])
    """
@typing.overload
def AddPatterns( sslib: SubstructLibrary, patterns: FPHolderBase, numThreads: int = 1) -> None:
    pass
def SubstructLibraryCanSerialize() -> bool:
    """
    SubstructLibraryCanSerialize() -> bool
        Returns True if the SubstructLibrary is serializable (requires boost serialization

        C++ signature :
            bool SubstructLibraryCanSerialize()
    """
