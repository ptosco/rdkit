"""Module containing functionality from the Avalon toolkit.

The functions currently exposed are:
  - GetCanonSmiles()   : return the canonical smiles for a molecule
  - GetAvalonFP()      : return the Avalon fingerprint for a molecule as
                         an RDKit ExplicitBitVector
  - GetAvalonCountFP()      : return the Avalon fingerprint for a molecule as
                              an RDKit SparseIntVector
  - Generate2DCoords() : use the Avalon coordinate generator to create
                         a set of 2D coordinates for a molecule
Each function can be called with either an RDKit molecule or some
molecule data as text (e.g. a SMILES or an MDL mol block).

See the individual docstrings for more information.
"""
from __future__ import annotations
import rdkit.Avalon.pyAvalonTools
import typing
import Boost.Python

__all__ = [
    "CheckMolecule",
    "CheckMoleculeString",
    "CloseCheckMolFiles",
    "Generate2DCoords",
    "GetAvalonCountFP",
    "GetAvalonFP",
    "GetAvalonFPAsWords",
    "GetCanonSmiles",
    "GetCheckMolLog",
    "InitializeCheckMol",
    "StruChkFlag",
    "StruChkResult",
    "avalonSSSBits",
    "avalonSimilarityBits"
]


class StruChkFlag(Boost.Python.enum, int):
    __slots__ = ()
    alias_conversion_failed = rdkit.Avalon.pyAvalonTools.StruChkFlag.alias_conversion_failed
    atom_check_failed = rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_check_failed
    atom_clash = rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_clash
    bad_molecule = rdkit.Avalon.pyAvalonTools.StruChkFlag.bad_molecule
    dubious_stereo_removed = rdkit.Avalon.pyAvalonTools.StruChkFlag.dubious_stereo_removed
    either_warning = rdkit.Avalon.pyAvalonTools.StruChkFlag.either_warning
    fragments_found = rdkit.Avalon.pyAvalonTools.StruChkFlag.fragments_found
    names = {'bad_molecule': rdkit.Avalon.pyAvalonTools.StruChkFlag.bad_molecule, 'alias_conversion_failed': rdkit.Avalon.pyAvalonTools.StruChkFlag.alias_conversion_failed, 'transformed': rdkit.Avalon.pyAvalonTools.StruChkFlag.transformed, 'fragments_found': rdkit.Avalon.pyAvalonTools.StruChkFlag.fragments_found, 'either_warning': rdkit.Avalon.pyAvalonTools.StruChkFlag.either_warning, 'stereo_error': rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_error, 'dubious_stereo_removed': rdkit.Avalon.pyAvalonTools.StruChkFlag.dubious_stereo_removed, 'atom_clash': rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_clash, 'atom_check_failed': rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_check_failed, 'size_check_failed': rdkit.Avalon.pyAvalonTools.StruChkFlag.size_check_failed, 'recharged': rdkit.Avalon.pyAvalonTools.StruChkFlag.recharged, 'stereo_forced_bad': rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_forced_bad, 'stereo_transformed': rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_transformed, 'template_transformed': rdkit.Avalon.pyAvalonTools.StruChkFlag.template_transformed}
    recharged = rdkit.Avalon.pyAvalonTools.StruChkFlag.recharged
    size_check_failed = rdkit.Avalon.pyAvalonTools.StruChkFlag.size_check_failed
    stereo_error = rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_error
    stereo_forced_bad = rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_forced_bad
    stereo_transformed = rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_transformed
    template_transformed = rdkit.Avalon.pyAvalonTools.StruChkFlag.template_transformed
    transformed = rdkit.Avalon.pyAvalonTools.StruChkFlag.transformed
    values = {1: rdkit.Avalon.pyAvalonTools.StruChkFlag.bad_molecule, 2: rdkit.Avalon.pyAvalonTools.StruChkFlag.alias_conversion_failed, 4: rdkit.Avalon.pyAvalonTools.StruChkFlag.transformed, 8: rdkit.Avalon.pyAvalonTools.StruChkFlag.fragments_found, 16: rdkit.Avalon.pyAvalonTools.StruChkFlag.either_warning, 32: rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_error, 64: rdkit.Avalon.pyAvalonTools.StruChkFlag.dubious_stereo_removed, 128: rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_clash, 256: rdkit.Avalon.pyAvalonTools.StruChkFlag.atom_check_failed, 512: rdkit.Avalon.pyAvalonTools.StruChkFlag.size_check_failed, 1024: rdkit.Avalon.pyAvalonTools.StruChkFlag.recharged, 2048: rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_forced_bad, 4096: rdkit.Avalon.pyAvalonTools.StruChkFlag.stereo_transformed, 8192: rdkit.Avalon.pyAvalonTools.StruChkFlag.template_transformed}
    pass
class StruChkResult(Boost.Python.enum, int):
    __slots__ = ()
    bad_set = rdkit.Avalon.pyAvalonTools.StruChkResult.bad_set
    names = {'success': rdkit.Avalon.pyAvalonTools.StruChkResult.success, 'bad_set': rdkit.Avalon.pyAvalonTools.StruChkResult.bad_set, 'transformed_set': rdkit.Avalon.pyAvalonTools.StruChkResult.transformed_set}
    success = rdkit.Avalon.pyAvalonTools.StruChkResult.success
    transformed_set = rdkit.Avalon.pyAvalonTools.StruChkResult.transformed_set
    values = {0: rdkit.Avalon.pyAvalonTools.StruChkResult.success, 2979: rdkit.Avalon.pyAvalonTools.StruChkResult.bad_set, 29788: rdkit.Avalon.pyAvalonTools.StruChkResult.transformed_set}
    pass
@typing.overload
def CheckMolecule( molstring: str, isSmiles: bool) -> tuple:
    """
    CheckMolecule( molstring: str, isSmiles: bool) -> tuple
        check a molecule passed in as a string.
        If the isSmiles argument is true, the string should represent the SMILES encoding
        of the molecule, otherwise it should be encoded as an MDL molfile.
        The first member of the return tuple contains the bit-encoded corrections made to the molecule.
        If possible, the molecule (corrected when appropriate) is returned as the second member of 
        the return tuple. Otherwise, None is returned.

        C++ signature :
            boost::python::tuple CheckMolecule(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool)

        C++ signature :
            boost::python::tuple CheckMolecule(RDKit::ROMol {lvalue})
    """
@typing.overload
def CheckMolecule( mol: object) -> tuple:
    pass
def CheckMoleculeString( molstring: str, isSmiles: bool) -> tuple:
    """
    CheckMoleculeString( molstring: str, isSmiles: bool) -> tuple
        check a molecule passed in as a string and returns the result as a string.
        If the isSmiles argument is true, the string should represent the SMILES encoding
        of the molecule, otherwise it should be encoded as an MDL molfile.
        The first member of the return tuple contains the bit-encoded corrections made to the molecule.
        If possible, a corrected CTAB for the molecule is returned as the second member of 
        the return tuple.

        C++ signature :
            boost::python::tuple CheckMoleculeString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool)
    """
def CloseCheckMolFiles() -> None:
    """
    CloseCheckMolFiles() -> None
        close open files used by molecule-checking functions.

        C++ signature :
            void CloseCheckMolFiles()
    """
@typing.overload
def Generate2DCoords( mol: object, clearConfs: bool = True) -> int:
    """
    Generate2DCoords( mol: object, clearConfs: bool = True) -> int
        Generates 2d coordinates for an RDKit molecule

        C++ signature :
            unsigned int Generate2DCoords(RDKit::ROMol {lvalue} [,bool=True])

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Generate2DCoords(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool)
    """
@typing.overload
def Generate2DCoords( molData: str, isSmiles: bool) -> str:
    pass
@typing.overload
def GetAvalonCountFP( mol: object, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407) -> object:
    """
    GetAvalonCountFP( mol: object, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407) -> object
        returns the Avalon count fingerprint for an RDKit molecule

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* GetAvalonCountFP(RDKit::ROMol [,unsigned int=512 [,bool=False [,unsigned int=15761407]]])

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* GetAvalonCountFP(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,unsigned int=15761407]]])
    """
@typing.overload
def GetAvalonCountFP( molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, bitFlags: int = 15761407) -> object:
    pass
@typing.overload
def GetAvalonFP( mol: object, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> object:
    """
    GetAvalonFP( mol: object, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> object
        returns the Avalon fingerprint for an RDKit molecule

        C++ signature :
            ExplicitBitVect* GetAvalonFP(RDKit::ROMol [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])

        C++ signature :
            ExplicitBitVect* GetAvalonFP(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])
    """
@typing.overload
def GetAvalonFP( molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> object:
    pass
@typing.overload
def GetAvalonFPAsWords( mol: object, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> list:
    """
    GetAvalonFPAsWords( mol: object, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> list
        returns the Avalon fingerprint for an RDKit molecule as a list of ints

        C++ signature :
            boost::python::list GetAvalonFPAsWords(RDKit::ROMol [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])

        C++ signature :
            boost::python::list GetAvalonFPAsWords(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,unsigned int=512 [,bool=False [,bool=False [,unsigned int=15761407]]]])
    """
@typing.overload
def GetAvalonFPAsWords( molData: str, isSmiles: bool, nBits: int = 512, isQuery: bool = False, resetVect: bool = False, bitFlags: int = 15761407) -> list:
    pass
@typing.overload
def GetCanonSmiles( mol: object, flags: int = -1) -> str:
    """
    GetCanonSmiles( mol: object, flags: int = -1) -> str
        returns canonical smiles for an RDKit molecule

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetCanonSmiles(RDKit::ROMol {lvalue} [,int=-1])

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetCanonSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,int=-1])
    """
@typing.overload
def GetCanonSmiles( molData: str, isSmiles: bool, flags: int = -1) -> str:
    pass
def GetCheckMolLog() -> str:
    """
    GetCheckMolLog() -> str
        Returns the Struchk log for the last molecules processed.

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetCheckMolLog()
    """
def InitializeCheckMol(options: str = '') -> int:
    """
    InitializeCheckMol(options: str = '') -> int
        initializes the structure checker.
        The argument should contain option lines separated by embedded newlines.An empty string will be used if the argument is omitted.An non-zero error code is returned in case of failure.

        C++ signature :
            int InitializeCheckMol([ std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=''])
    """
avalonSSSBits = 32767
avalonSimilarityBits = 15761407
