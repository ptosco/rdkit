from __future__ import annotations
import rdkit.Chem.rdMHFPFingerprint
import typing
import Boost.Python

__all__ = [
    "MHFPEncoder"
]


class MHFPEncoder(Boost.Python.instance):
    @staticmethod
    def CreateShinglingFromMol( arg1: MHFPEncoder, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        CreateShinglingFromMol( arg1: MHFPEncoder, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > CreateShinglingFromMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    @staticmethod
    def CreateShinglingFromSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        CreateShinglingFromSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > CreateShinglingFromSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    @staticmethod
    def Distance( arg1: _vectj, arg2: _vectj) -> float: 
        """
        Distance( arg1: _vectj, arg2: _vectj) -> float

            C++ signature :
                double Distance(std::vector<unsigned int, std::allocator<unsigned int> >,std::vector<unsigned int, std::allocator<unsigned int> >)
        """
    @staticmethod
    def EncodeMol( arg1: MHFPEncoder, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectj: 
        """
        EncodeMol( arg1: MHFPEncoder, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectj
            Creates a MHFP vector from an RDKit Mol instance.

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > EncodeMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    @staticmethod
    def EncodeMolsBulk( arg1: MHFPEncoder, mols: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectSt6vectorIjSaIjEE: 
        """
        EncodeMolsBulk( arg1: MHFPEncoder, mols: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectSt6vectorIjSaIjEE
            Creates a MHFP vector from a list of RDKit Mol instances.

            C++ signature :
                std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > > EncodeMolsBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    @staticmethod
    def EncodeSECFPMol( arg1: MHFPEncoder, smiles: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect: 
        """
        EncodeSECFPMol( arg1: MHFPEncoder, smiles: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect
            Creates a SECFP binary vector from an RDKit Mol instance.

            C++ signature :
                ExplicitBitVect EncodeSECFPMol(RDKit::MHFPFingerprints::MHFPEncoder*,RDKit::ROMol [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
    @staticmethod
    def EncodeSECFPMolsBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> object: 
        """
        EncodeSECFPMolsBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> object
            Creates a SECFP binary vector from a list of RDKit Mol instances.

            C++ signature :
                std::vector<ExplicitBitVect, std::allocator<ExplicitBitVect> > EncodeSECFPMolsBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
    @staticmethod
    def EncodeSECFPSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect: 
        """
        EncodeSECFPSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect
            Creates a SECFP binary vector from a SMILES string.

            C++ signature :
                ExplicitBitVect EncodeSECFPSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
    @staticmethod
    def EncodeSECFPSmilesBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> object: 
        """
        EncodeSECFPSmilesBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> object
            Creates a SECFP binary vector from a list of SMILES strings.

            C++ signature :
                std::vector<ExplicitBitVect, std::allocator<ExplicitBitVect> > EncodeSECFPSmilesBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1 [,unsigned long=2048]]]]]])
        """
    @staticmethod
    def EncodeSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectj: 
        """
        EncodeSmiles( arg1: MHFPEncoder, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectj
            Creates a MHFP vector from a SMILES string.

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > EncodeSmiles(RDKit::MHFPFingerprints::MHFPEncoder*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    @staticmethod
    def EncodeSmilesBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectSt6vectorIjSaIjEE: 
        """
        EncodeSmilesBulk( arg1: MHFPEncoder, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectSt6vectorIjSaIjEE
            Creates a MHFP vector from a list of SMILES strings.

            C++ signature :
                std::vector<std::vector<unsigned int, std::allocator<unsigned int> >, std::allocator<std::vector<unsigned int, std::allocator<unsigned int> > > > EncodeSmilesBulk(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue} [,unsigned char=3 [,bool=True [,bool=False [,bool=False [,unsigned char=1]]]]])
        """
    def FromArray(self, vec: list) -> _vectj: 
        """
        FromArray( self: MHFPEncoder, vec: list) -> _vectj
            Creates a MHFP vector from a list of unsigned integers.

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > FromArray(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue})
        """
    def FromStringArray(self, vec: list) -> _vectj: 
        """
        FromStringArray( self: MHFPEncoder, vec: list) -> _vectj
            Creates a MHFP vector from a list of arbitrary strings.

            C++ signature :
                std::vector<unsigned int, std::allocator<unsigned int> > FromStringArray(RDKit::MHFPFingerprints::MHFPEncoder*,boost::python::list {lvalue})
        """
    @staticmethod
    def __init__( arg1: object, arg2: int, arg3: int) -> None: 
        """
        __init__( arg1: object, arg2: int, arg3: int) -> None

            C++ signature :
                void __init__(_object* [,unsigned int [,unsigned int]])
        """
    __instance_size__ = 96
    pass
