"""Module containing RDKit functionality for working with molecular file formats."""
from __future__ import annotations
import rdkit.Chem.rdmolfiles
import typing
import Boost.Python

__all__ = [
    "AddMetadataToPNGFile",
    "AddMetadataToPNGString",
    "AtomFromSmarts",
    "AtomFromSmiles",
    "BondFromSmarts",
    "BondFromSmiles",
    "CXSmilesFields",
    "CanonicalRankAtoms",
    "CanonicalRankAtomsInFragment",
    "CanonicalizeEnhancedStereo",
    "CreateAtomBoolPropertyList",
    "CreateAtomDoublePropertyList",
    "CreateAtomIntPropertyList",
    "CreateAtomStringPropertyList",
    "ForwardSDMolSupplier",
    "MaeMolSupplier",
    "MaeWriter",
    "MetadataFromPNGFile",
    "MetadataFromPNGString",
    "MolFragmentToCXSmarts",
    "MolFragmentToCXSmiles",
    "MolFragmentToSmarts",
    "MolFragmentToSmiles",
    "MolFromFASTA",
    "MolFromHELM",
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
    "MolToCMLBlock",
    "MolToCMLFile",
    "MolToCXSmarts",
    "MolToCXSmiles",
    "MolToFASTA",
    "MolToHELM",
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
    "MultithreadedSDMolSupplier",
    "MultithreadedSmilesMolSupplier",
    "PDBWriter",
    "SDMolSupplier",
    "SDWriter",
    "SmartsParserParams",
    "SmilesMolSupplier",
    "SmilesMolSupplierFromText",
    "SmilesParserParams",
    "SmilesWriteParams",
    "SmilesWriter",
    "TDTMolSupplier",
    "TDTWriter"
]


class CXSmilesFields(Boost.Python.enum, int):
    CX_ALL = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL
    CX_ALL_BUT_COORDS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL_BUT_COORDS
    CX_ATOM_LABELS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_LABELS
    CX_ATOM_PROPS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_PROPS
    CX_BOND_CFG = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_BOND_CFG
    CX_COORDS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_COORDS
    CX_ENHANCEDSTEREO = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ENHANCEDSTEREO
    CX_LINKNODES = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_LINKNODES
    CX_MOLFILE_VALUES = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_MOLFILE_VALUES
    CX_NONE = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_NONE
    CX_POLYMER = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_POLYMER
    CX_RADICALS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_RADICALS
    CX_SGROUPS = rdkit.Chem.rdmolfiles.CXSmilesFields.CX_SGROUPS
    __slots__ = ()
    names = {'CX_NONE': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_NONE, 'CX_ATOM_LABELS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_LABELS, 'CX_MOLFILE_VALUES': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_MOLFILE_VALUES, 'CX_COORDS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_COORDS, 'CX_RADICALS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_RADICALS, 'CX_ATOM_PROPS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_PROPS, 'CX_LINKNODES': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_LINKNODES, 'CX_ENHANCEDSTEREO': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ENHANCEDSTEREO, 'CX_SGROUPS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_SGROUPS, 'CX_POLYMER': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_POLYMER, 'CX_BOND_CFG': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_BOND_CFG, 'CX_ALL': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL, 'CX_ALL_BUT_COORDS': rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL_BUT_COORDS}
    values = {0: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_NONE, 1: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_LABELS, 2: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_MOLFILE_VALUES, 4: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_COORDS, 8: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_RADICALS, 16: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ATOM_PROPS, 32: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_LINKNODES, 64: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ENHANCEDSTEREO, 128: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_SGROUPS, 256: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_POLYMER, 512: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_BOND_CFG, 2147483647: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL, 2147483643: rdkit.Chem.rdmolfiles.CXSmilesFields.CX_ALL_BUT_COORDS}
    pass
class ForwardSDMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from file-like object containing SD data.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = ForwardSDMolSupplier(file('in.sdf'))
           >>> for mol in suppl:
           ...    if mol is not None: mol.GetNumAtoms()

        2) we can also read from compressed files: 

           >>> import gzip
           >>> suppl = ForwardSDMolSupplier(gzip.open('in.sdf.gz'))
           >>> for mol in suppl:
           ...   if mol is not None: print mol.GetNumAtoms()

      Properties in the SD file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """
    @staticmethod
    def GetEOFHitOnRead( arg1: ForwardSDMolSupplier) -> bool: 
        """
        GetEOFHitOnRead( arg1: ForwardSDMolSupplier) -> bool
            Returns whether or EOF was hit while parsing the previous entry.
            

            C++ signature :
                bool GetEOFHitOnRead((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
    @staticmethod
    def GetProcessPropertyLists( arg1: ForwardSDMolSupplier) -> bool: 
        """
        GetProcessPropertyLists( arg1: ForwardSDMolSupplier) -> bool
            returns whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                bool GetProcessPropertyLists((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
    @staticmethod
    def SetProcessPropertyLists( arg1: ForwardSDMolSupplier, arg2: bool) -> None: 
        """
        SetProcessPropertyLists( arg1: ForwardSDMolSupplier, arg2: bool) -> None
            sets whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                void SetProcessPropertyLists((anonymous namespace)::LocalForwardSDMolSupplier {lvalue},bool)
        """
    @staticmethod
    def __enter__( arg1: ForwardSDMolSupplier) -> ForwardSDMolSupplier: 
        """
        __enter__( arg1: ForwardSDMolSupplier) -> ForwardSDMolSupplier

            C++ signature :
                (anonymous namespace)::LocalForwardSDMolSupplier* __enter__((anonymous namespace)::LocalForwardSDMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: ForwardSDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: ForwardSDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__((anonymous namespace)::LocalForwardSDMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileobj: object, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: 
        """
        __init__( arg1: object, fileobj: object, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None

            C++ signature :
                void __init__(_object*,boost::python::api::object {lvalue} [,bool=True [,bool=True [,bool=True]]])

            C++ signature :
                void __init__(_object*,boost_adaptbx::python::streambuf {lvalue} [,bool=True [,bool=True [,bool=True]]])

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, streambuf: streambuf, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, filename: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...
    @staticmethod
    def __iter__( arg1: ForwardSDMolSupplier) -> ForwardSDMolSupplier: 
        """
        __iter__( arg1: ForwardSDMolSupplier) -> ForwardSDMolSupplier

            C++ signature :
                (anonymous namespace)::LocalForwardSDMolSupplier* __iter__((anonymous namespace)::LocalForwardSDMolSupplier*)
        """
    @staticmethod
    def __next__( arg1: ForwardSDMolSupplier) -> Mol: 
        """
        __next__( arg1: ForwardSDMolSupplier) -> Mol
            Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__((anonymous namespace)::LocalForwardSDMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: ForwardSDMolSupplier) -> bool: 
        """
        atEnd( arg1: ForwardSDMolSupplier) -> bool
            Returns whether or not we have hit EOF.
            

            C++ signature :
                bool atEnd((anonymous namespace)::LocalForwardSDMolSupplier {lvalue})
        """
    pass
class MaeMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from file-like object containing Maestro data.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = MaeMolSupplier(file('in.mae'))
           >>> for mol in suppl:
           ...    if mol is not None: mol.GetNumAtoms()

        2) we can also read from compressed files: 

           >>> import gzip
           >>> suppl = MaeMolSupplier(gzip.open('in.maegz'))
           >>> for mol in suppl:
           ...   if mol is not None: print mol.GetNumAtoms()

      Properties in the Maestro file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """
    @staticmethod
    def SetData( arg1: MaeMolSupplier, data: str, sanitize: bool = True, removeHs: bool = True) -> None: 
        """
        SetData( arg1: MaeMolSupplier, data: str, sanitize: bool = True, removeHs: bool = True) -> None
            Sets the text to be parsed

            C++ signature :
                void SetData((anonymous namespace)::LocalMaeMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
        """
    @staticmethod
    def __enter__( arg1: MaeMolSupplier) -> MaeMolSupplier: 
        """
        __enter__( arg1: MaeMolSupplier) -> MaeMolSupplier

            C++ signature :
                (anonymous namespace)::LocalMaeMolSupplier* __enter__((anonymous namespace)::LocalMaeMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: MaeMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: MaeMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__((anonymous namespace)::LocalMaeMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def __getitem__( arg1: MaeMolSupplier, arg2: int) -> Mol: 
        """
        __getitem__( arg1: MaeMolSupplier, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__((anonymous namespace)::LocalMaeMolSupplier*,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,boost::python::api::object {lvalue} [,bool=True [,bool=True]])

            C++ signature :
                void __init__(_object*,boost_adaptbx::python::streambuf {lvalue} [,bool=True [,bool=True]])

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileobj: object, sanitize: bool = True, removeHs: bool = True) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, streambuf: streambuf, sanitize: bool = True, removeHs: bool = True) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, filename: str, sanitize: bool = True, removeHs: bool = True) -> None: ...
    @staticmethod
    def __iter__( arg1: MaeMolSupplier) -> MaeMolSupplier: 
        """
        __iter__( arg1: MaeMolSupplier) -> MaeMolSupplier

            C++ signature :
                (anonymous namespace)::LocalMaeMolSupplier* __iter__((anonymous namespace)::LocalMaeMolSupplier*)
        """
    @staticmethod
    def __len__( arg1: MaeMolSupplier) -> int: 
        """
        __len__( arg1: MaeMolSupplier) -> int

            C++ signature :
                unsigned int __len__((anonymous namespace)::LocalMaeMolSupplier {lvalue})
        """
    @staticmethod
    def __next__( arg1: MaeMolSupplier) -> Mol: 
        """
        __next__( arg1: MaeMolSupplier) -> Mol
            Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__((anonymous namespace)::LocalMaeMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: MaeMolSupplier) -> bool: 
        """
        atEnd( arg1: MaeMolSupplier) -> bool
            Returns whether or not we have hit EOF.
            

            C++ signature :
                bool atEnd((anonymous namespace)::LocalMaeMolSupplier {lvalue})
        """
    @staticmethod
    def reset( arg1: MaeMolSupplier) -> None: 
        """
        reset( arg1: MaeMolSupplier) -> None
            Resets our position in the file to the beginning.
            

            C++ signature :
                void reset((anonymous namespace)::LocalMaeMolSupplier {lvalue})
        """
    __instance_size__ = 144
    pass
class MaeWriter(Boost.Python.instance):
    """
    An experimental class for writing molecules to Maestro files.

      Usage examples:

        1) writing to a named file:

           >>> writer = MaeWriter('out.mae')
           >>> for mol in list_of_mols:
           ...    writer.write(mol)

        2) writing to a file-like object: 

           >>> import gzip
           >>> outf=gzip.open('out.mae.gz','wt+')
           >>> writer = MaeWriter(outf)
           >>> for mol in list_of_mols:
           ...   writer.write(mol)
           >>> writer.close()
           >>> outf.close()

      By default all non-private molecule, atom and bond properties are written
      to the Maestro file. This can be changed using the SetProps method:

           >>> writer = MaeWriter('out.mae')
           >>> writer.SetProps(['prop1','prop2'])

      Properties that are specified, but are not present will be ignored.

      Kekulization is mandatory, as the Maestro format does not have
      the concept of an aromatic bond

      As this is an experimental writer, many features are not supported yet,
      e.g. chirality and bond stereo labels, stereo groups, substance groups,
      isotopes, or even dummy atoms. Note that these aren't supported by
      MaeMolSupplier either.

     
    """
    @staticmethod
    def GetText( mol: Mol, confId: int = -1, props_list: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE = _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE()) -> str: 
        """
        GetText( mol: Mol, confId: int = -1, props_list: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE = _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE()) -> str
            returns the Maestro ct block text for a molecule

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetText(RDKit::ROMol [,int=-1 [,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >=<rdkit.rdBase._vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE object at 0x7f2c65395dc0>]])
        """
    @staticmethod
    def NumMols( arg1: MaeWriter) -> int: 
        """
        NumMols( arg1: MaeWriter) -> int
            Returns the number of molecules written so far.
            
            

            C++ signature :
                unsigned int NumMols(RDKit::LocalMaeWriter {lvalue})
        """
    def SetProps(self, props_list: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE) -> None: 
        """
        SetProps( self: MaeWriter, props_list: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE) -> None
            Sets the atom and mol properties to be written to the output file
            
              ARGUMENTS:
            
                - props: a list or tuple of atom and mol property names
            
            

            C++ signature :
                void SetProps(RDKit::LocalMaeWriter {lvalue},std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)
        """
    @staticmethod
    def __enter__( arg1: MaeWriter) -> MaeWriter: 
        """
        __enter__( arg1: MaeWriter) -> MaeWriter

            C++ signature :
                RDKit::LocalMaeWriter* __enter__(RDKit::LocalMaeWriter*)
        """
    @staticmethod
    def __exit__( arg1: MaeWriter, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: MaeWriter, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::LocalMaeWriter*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileobj: object) -> None: 
        """
        __init__( arg1: object, fileobj: object) -> None

            C++ signature :
                void __init__(_object*,boost::python::api::object {lvalue})

            C++ signature :
                void __init__(_object*,boost_adaptbx::python::streambuf {lvalue})

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, streambuf: streambuf) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, filename: str) -> None: ...
    @staticmethod
    def close( arg1: MaeWriter) -> None: 
        """
        close( arg1: MaeWriter) -> None
            Flushes the output file and closes it. The Writer cannot be used after this.
            
            

            C++ signature :
                void close(RDKit::LocalMaeWriter {lvalue})
        """
    @staticmethod
    def flush( arg1: MaeWriter) -> None: 
        """
        flush( arg1: MaeWriter) -> None
            Flushes the output file (forces the disk file to be updated).
            
            

            C++ signature :
                void flush(RDKit::LocalMaeWriter {lvalue})
        """
    def write(self, mol: Mol, confId: int = -1) -> None: 
        """
        write( self: MaeWriter, mol: Mol, confId: int = -1) -> None
            Writes a molecule to the output file.
            
              ARGUMENTS:
            
                - mol: the Mol to be written
                - confId: (optional) ID of the conformation to write
            
            

            C++ signature :
                void write(RDKit::LocalMaeWriter {lvalue},RDKit::ROMol [,int=-1])
        """
    pass
class MultithreadedSDMolSupplier(Boost.Python.instance):
    """
    A class which concurrently supplies molecules from a text file.
      Please note that this class is still a bit experimental and the API may
      change in future releases.

      Usage examples:

        1) Lazy evaluation: the molecules might not be constructed until we ask for them:

           >>> suppl = MultithreadedSDMolSupplier('in.sdf')
           >>> for mol in suppl:
           ...    if(mol):
           ...      mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = MultithreadedSDMolSupplier('in.sdf')
           >>> while (!suppl.atEnd()):
           ...    mol = next(mol)
           ...    if(mol):
           ...      mol.GetNumAtoms()
    """
    @staticmethod
    def GetLastItemText( arg1: MultithreadedSDMolSupplier) -> str: 
        """
        GetLastItemText( arg1: MultithreadedSDMolSupplier) -> str
            Returns the text for the last extracted item.
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetLastItemText(RDKit::MultithreadedSDMolSupplier*)
        """
    @staticmethod
    def GetLastRecordId( arg1: MultithreadedSDMolSupplier) -> int: 
        """
        GetLastRecordId( arg1: MultithreadedSDMolSupplier) -> int
            Returns the record id for the last extracted item.
            

            C++ signature :
                unsigned int GetLastRecordId(RDKit::MultithreadedSDMolSupplier*)
        """
    @staticmethod
    def GetProcessPropertyLists( arg1: MultithreadedSDMolSupplier) -> bool: 
        """
        GetProcessPropertyLists( arg1: MultithreadedSDMolSupplier) -> bool
            returns whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                bool GetProcessPropertyLists(RDKit::MultithreadedSDMolSupplier {lvalue})
        """
    @staticmethod
    def SetProcessPropertyLists( arg1: MultithreadedSDMolSupplier, arg2: bool) -> None: 
        """
        SetProcessPropertyLists( arg1: MultithreadedSDMolSupplier, arg2: bool) -> None
            sets whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                void SetProcessPropertyLists(RDKit::MultithreadedSDMolSupplier {lvalue},bool)
        """
    @staticmethod
    def __enter__( arg1: MultithreadedSDMolSupplier) -> MultithreadedSDMolSupplier: 
        """
        __enter__( arg1: MultithreadedSDMolSupplier) -> MultithreadedSDMolSupplier

            C++ signature :
                RDKit::MultithreadedSDMolSupplier* __enter__(RDKit::MultithreadedSDMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: MultithreadedSDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: MultithreadedSDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::MultithreadedSDMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True [,unsigned int=1 [,unsigned long=5 [,unsigned long=5]]]]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True, numWriterThreads: int = 1, sizeInputQueue: int = 5, sizeOutputQueue: int = 5) -> None: ...
    @staticmethod
    def __iter__( arg1: MultithreadedSDMolSupplier) -> MultithreadedSDMolSupplier: 
        """
        __iter__( arg1: MultithreadedSDMolSupplier) -> MultithreadedSDMolSupplier

            C++ signature :
                RDKit::MultithreadedSDMolSupplier* __iter__(RDKit::MultithreadedSDMolSupplier*)
        """
    @staticmethod
    def __next__( arg1: MultithreadedSDMolSupplier) -> Mol: 
        """
        __next__( arg1: MultithreadedSDMolSupplier) -> Mol
            Returns the next molecule in the file. Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::MultithreadedSDMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: MultithreadedSDMolSupplier) -> bool: 
        """
        atEnd( arg1: MultithreadedSDMolSupplier) -> bool
            Returns true if we have read all records else false.
            

            C++ signature :
                bool atEnd(RDKit::MultithreadedSDMolSupplier {lvalue})
        """
    __instance_size__ = 184
    pass
class MultithreadedSmilesMolSupplier(Boost.Python.instance):
    """
    A class which concurrently supplies molecules from a text file.
      Please note that this class is still a bit experimental and the API may
      change in future releases.

      Usage examples:

        1) Lazy evaluation: the molecules might not be constructed until we ask for them:

           >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    if(mol):
           ...      mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = MultithreadedSmilesMolSupplier('in.smi')
           >>> while (!suppl.atEnd()):
           ...    mol = next(mol)
           ...    if(mol):
           ...      mol.GetNumAtoms()
    """
    @staticmethod
    def GetLastItemText( arg1: MultithreadedSmilesMolSupplier) -> str: 
        """
        GetLastItemText( arg1: MultithreadedSmilesMolSupplier) -> str
            Returns the text for the last extracted item.
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetLastItemText(RDKit::MultithreadedSmilesMolSupplier*)
        """
    @staticmethod
    def GetLastRecordId( arg1: MultithreadedSmilesMolSupplier) -> int: 
        """
        GetLastRecordId( arg1: MultithreadedSmilesMolSupplier) -> int
            Returns the record id for the last extracted item.
            

            C++ signature :
                unsigned int GetLastRecordId(RDKit::MultithreadedSmilesMolSupplier*)
        """
    @staticmethod
    def __enter__( arg1: MultithreadedSmilesMolSupplier) -> MultithreadedSmilesMolSupplier: 
        """
        __enter__( arg1: MultithreadedSmilesMolSupplier) -> MultithreadedSmilesMolSupplier

            C++ signature :
                RDKit::MultithreadedSmilesMolSupplier* __enter__(RDKit::MultithreadedSmilesMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: MultithreadedSmilesMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: MultithreadedSmilesMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::MultithreadedSmilesMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' \t' [,int=0 [,int=1 [,bool=True [,bool=True [,unsigned int=1 [,unsigned long=5 [,unsigned long=5]]]]]]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, delimiter: str = '\t', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True, numWriterThreads: int = 1, sizeInputQueue: int = 5, sizeOutputQueue: int = 5) -> None: ...
    @staticmethod
    def __iter__( arg1: MultithreadedSmilesMolSupplier) -> MultithreadedSmilesMolSupplier: 
        """
        __iter__( arg1: MultithreadedSmilesMolSupplier) -> MultithreadedSmilesMolSupplier

            C++ signature :
                RDKit::MultithreadedSmilesMolSupplier* __iter__(RDKit::MultithreadedSmilesMolSupplier*)
        """
    @staticmethod
    def __next__( arg1: MultithreadedSmilesMolSupplier) -> Mol: 
        """
        __next__( arg1: MultithreadedSmilesMolSupplier) -> Mol
            Returns the next molecule in the file. Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::MultithreadedSmilesMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: MultithreadedSmilesMolSupplier) -> bool: 
        """
        atEnd( arg1: MultithreadedSmilesMolSupplier) -> bool
            Returns true if we have read all records else false.
            

            C++ signature :
                bool atEnd(RDKit::MultithreadedSmilesMolSupplier {lvalue})
        """
    __instance_size__ = 248
    pass
class PDBWriter(Boost.Python.instance):
    """
    A class for writing molecules to PDB files.
    """
    @staticmethod
    def NumMols( arg1: PDBWriter) -> int: 
        """
        NumMols( arg1: PDBWriter) -> int
            Returns the number of molecules written so far.
            
            

            C++ signature :
                unsigned int NumMols(RDKit::PDBWriter {lvalue})
        """
    @staticmethod
    def __enter__( arg1: PDBWriter) -> PDBWriter: 
        """
        __enter__( arg1: PDBWriter) -> PDBWriter

            C++ signature :
                RDKit::PDBWriter* __enter__(RDKit::PDBWriter*)
        """
    @staticmethod
    def __exit__( arg1: PDBWriter, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: PDBWriter, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::PDBWriter*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileObj: object, flavor: int = 0) -> object: 
        """
        __init__( arg1: object, fileObj: object, flavor: int = 0) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object {lvalue} [,unsigned int=0])

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,unsigned int=0])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, flavor: int = 0) -> None: ...
    @staticmethod
    def close( arg1: PDBWriter) -> None: 
        """
        close( arg1: PDBWriter) -> None
            Flushes the output file and closes it. The Writer cannot be used after this.
            
            

            C++ signature :
                void close(RDKit::PDBWriter {lvalue})
        """
    @staticmethod
    def flush( arg1: PDBWriter) -> None: 
        """
        flush( arg1: PDBWriter) -> None
            Flushes the output file (forces the disk file to be updated).
            
            

            C++ signature :
                void flush(RDKit::PDBWriter {lvalue})
        """
    def write(self, mol: Mol, confId: int = -1) -> None: 
        """
        write( self: PDBWriter, mol: Mol, confId: int = -1) -> None
            Writes a molecule to the output file.
            
              ARGUMENTS:
            
                - mol: the Mol to be written
                - confId: (optional) ignored 
            
            

            C++ signature :
                void write(RDKit::PDBWriter {lvalue},RDKit::ROMol [,int=-1])
        """
    pass
class SDMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from an SD file.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = SDMolSupplier('in.sdf')
           >>> for mol in suppl:
           ...    mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = SDMolSupplier('in.sdf')
           >>> mol1 = next(suppl)
           >>> mol2 = next(suppl)
           >>> suppl.reset()
           >>> mol3 = next(suppl)
           # mol3 and mol1 are the same:
           >>> MolToSmiles(mol3)==MolToSmiles(mol1)

        3) Random Access:

           >>> suppl = SDMolSupplier('in.sdf')
           >>> mol1 = suppl[0] 
           >>> mol2 = suppl[1] 
           # NOTE: this will generate an IndexError if the supplier doesn't have that many
           molecules.

        4) Random Access 2:  looping over all molecules 

           >>> suppl = SDMolSupplier('in.sdf')
           >>> nMols = len(suppl)
           >>> for i in range(nMols):
           ...   suppl[i].GetNumAtoms()

      Properties in the SD file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """
    def GetItemText(self, index: int) -> str: 
        """
        GetItemText( self: SDMolSupplier, index: int) -> str
            returns the text for an item

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::SDMolSupplier {lvalue},unsigned int)
        """
    @staticmethod
    def GetProcessPropertyLists( arg1: SDMolSupplier) -> bool: 
        """
        GetProcessPropertyLists( arg1: SDMolSupplier) -> bool
            returns whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                bool GetProcessPropertyLists(RDKit::SDMolSupplier {lvalue})
        """
    def SetData(self, data: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: 
        """
        SetData( self: SDMolSupplier, data: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None
            Sets the text to be parsed

            C++ signature :
                void SetData(RDKit::SDMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    def SetProcessPropertyLists( arg1: SDMolSupplier, arg2: bool) -> None: 
        """
        SetProcessPropertyLists( arg1: SDMolSupplier, arg2: bool) -> None
            sets whether or not any property lists that are present will be processed when reading molecules

            C++ signature :
                void SetProcessPropertyLists(RDKit::SDMolSupplier {lvalue},bool)
        """
    def _SetStreamIndices(self, locs: object) -> None: 
        """
        _SetStreamIndices( self: SDMolSupplier, locs: object) -> None
            Sets the locations of mol beginnings in the input stream. Be *very* careful with this method.

            C++ signature :
                void _SetStreamIndices(RDKit::SDMolSupplier {lvalue},boost::python::api::object)
        """
    @staticmethod
    def __enter__( arg1: SDMolSupplier) -> SDMolSupplier: 
        """
        __enter__( arg1: SDMolSupplier) -> SDMolSupplier

            C++ signature :
                RDKit::SDMolSupplier* __enter__(RDKit::SDMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: SDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: SDMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::SDMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def __getitem__( arg1: SDMolSupplier, arg2: int) -> Mol: 
        """
        __getitem__( arg1: SDMolSupplier, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__(RDKit::SDMolSupplier*,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True [,bool=True]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, sanitize: bool = True, removeHs: bool = True, strictParsing: bool = True) -> None: ...
    @staticmethod
    def __iter__( arg1: SDMolSupplier) -> SDMolSupplier: 
        """
        __iter__( arg1: SDMolSupplier) -> SDMolSupplier

            C++ signature :
                RDKit::SDMolSupplier* __iter__(RDKit::SDMolSupplier*)
        """
    @staticmethod
    def __len__( arg1: SDMolSupplier) -> int: 
        """
        __len__( arg1: SDMolSupplier) -> int

            C++ signature :
                unsigned int __len__(RDKit::SDMolSupplier {lvalue})
        """
    @staticmethod
    def __next__( arg1: SDMolSupplier) -> Mol: 
        """
        __next__( arg1: SDMolSupplier) -> Mol
            Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::SDMolSupplier*)
        """
    @staticmethod
    def atEnd( arg1: SDMolSupplier) -> bool: 
        """
        atEnd( arg1: SDMolSupplier) -> bool
            Returns whether or not we have hit EOF.
            

            C++ signature :
                bool atEnd(RDKit::SDMolSupplier {lvalue})
        """
    @staticmethod
    def reset( arg1: SDMolSupplier) -> None: 
        """
        reset( arg1: SDMolSupplier) -> None
            Resets our position in the file to the beginning.
            

            C++ signature :
                void reset(RDKit::SDMolSupplier {lvalue})
        """
    __instance_size__ = 88
    pass
class SDWriter(Boost.Python.instance):
    """
    A class for writing molecules to SD files.

      Usage examples:

        1) writing to a named file:

           >>> writer = SDWriter('out.sdf')
           >>> for mol in list_of_mols:
           ...    writer.write(mol)

        2) writing to a file-like object: 

           >>> import gzip
           >>> outf=gzip.open('out.sdf.gz','wt+')
           >>> writer = SDWriter(outf)
           >>> for mol in list_of_mols:
           ...   writer.write(mol)
           >>> writer.close()
           >>> outf.close()

      By default all non-private molecular properties are written to the SD file.
      This can be changed using the SetProps method:

           >>> writer = SDWriter('out.sdf')
           >>> writer.SetProps(['prop1','prop2'])
    """
    @staticmethod
    def GetForceV3000( arg1: SDWriter) -> bool: 
        """
        GetForceV3000( arg1: SDWriter) -> bool
            Returns whether or not V3000 mol file writing is being forced.
            
            

            C++ signature :
                bool GetForceV3000(RDKit::SDWriter {lvalue})
        """
    @staticmethod
    def GetKekulize( arg1: SDWriter) -> bool: 
        """
        GetKekulize( arg1: SDWriter) -> bool
            Returns whether or not molecules are kekulized on writing.
            
            

            C++ signature :
                bool GetKekulize(RDKit::SDWriter {lvalue})
        """
    @staticmethod
    def GetText( mol: Mol, confId: int = -1, kekulize: bool = True, force_v3000: bool = False, molid: int = -1) -> str: 
        """
        GetText( mol: Mol, confId: int = -1, kekulize: bool = True, force_v3000: bool = False, molid: int = -1) -> str
            returns the SD text for a molecule

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetText(RDKit::ROMol [,int=-1 [,bool=True [,bool=False [,int=-1]]]])
        """
    @staticmethod
    def NumMols( arg1: SDWriter) -> int: 
        """
        NumMols( arg1: SDWriter) -> int
            Returns the number of molecules written so far.
            
            

            C++ signature :
                unsigned int NumMols(RDKit::SDWriter {lvalue})
        """
    @staticmethod
    def SetForceV3000( arg1: SDWriter, arg2: bool) -> None: 
        """
        SetForceV3000( arg1: SDWriter, arg2: bool) -> None
            Sets whether or not V3000 mol file writing is being forced.
            
            

            C++ signature :
                void SetForceV3000(RDKit::SDWriter {lvalue},bool)
        """
    @staticmethod
    def SetKekulize( arg1: SDWriter, arg2: bool) -> None: 
        """
        SetKekulize( arg1: SDWriter, arg2: bool) -> None
            Sets whether or not molecules are kekulized on writing.
            
            

            C++ signature :
                void SetKekulize(RDKit::SDWriter {lvalue},bool)
        """
    @staticmethod
    def SetProps( arg1: SDWriter, arg2: object) -> None: 
        """
        SetProps( arg1: SDWriter, arg2: object) -> None
            Sets the properties to be written to the output file
            
              ARGUMENTS:
            
                - props: a list or tuple of property names
            
            

            C++ signature :
                void SetProps(RDKit::SDWriter {lvalue},boost::python::api::object)
        """
    @staticmethod
    def __enter__( arg1: SDWriter) -> SDWriter: 
        """
        __enter__( arg1: SDWriter) -> SDWriter

            C++ signature :
                RDKit::SDWriter* __enter__(RDKit::SDWriter*)
        """
    @staticmethod
    def __exit__( arg1: SDWriter, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: SDWriter, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::SDWriter*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: object) -> object: 
        """
        __init__( arg1: object, arg2: object) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object {lvalue})

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str) -> None: ...
    @staticmethod
    def close( arg1: SDWriter) -> None: 
        """
        close( arg1: SDWriter) -> None
            Flushes the output file and closes it. The Writer cannot be used after this.
            
            

            C++ signature :
                void close(RDKit::SDWriter {lvalue})
        """
    @staticmethod
    def flush( arg1: SDWriter) -> None: 
        """
        flush( arg1: SDWriter) -> None
            Flushes the output file (forces the disk file to be updated).
            
            

            C++ signature :
                void flush(RDKit::SDWriter {lvalue})
        """
    def write(self, mol: Mol, confId: int = -1) -> None: 
        """
        write( self: SDWriter, mol: Mol, confId: int = -1) -> None
            Writes a molecule to the output file.
            
              ARGUMENTS:
            
                - mol: the Mol to be written
                - confId: (optional) ID of the conformation to write
            
            

            C++ signature :
                void write(RDKit::SDWriter {lvalue},RDKit::ROMol {lvalue} [,int=-1])
        """
    pass
class SmartsParserParams(Boost.Python.instance):
    """
    Parameters controlling SMARTS Parsing
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def allowCXSMILES(self) -> None:
        """
        controls whether or not the CXSMILES extensions are parsed

        :type: None
        """
    @property
    def debugParse(self) -> None:
        """
        controls the amount of debugging information produced

        :type: None
        """
    @property
    def mergeHs(self) -> None:
        """
        toggles merging H atoms in the SMARTS into neighboring atoms

        :type: None
        """
    @property
    def parseName(self) -> None:
        """
        controls whether or not the molecule name is also parsed

        :type: None
        """
    @property
    def strictCXSMILES(self) -> None:
        """
        controls whether or not problems in CXSMILES parsing causes molecule parsing to fail

        :type: None
        """
    __instance_size__ = 48
    pass
class SmilesMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from a text file.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> mol1 = next(suppl)
           >>> mol2 = next(suppl)
           >>> suppl.reset()
           >>> mol3 = next(suppl)
           # mol3 and mol1 are the same:
           >>> MolToSmiles(mol3)==MolToSmiles(mol1)

        3) Random Access:  all molecules are constructed as soon as we ask for the
           length:

           >>> suppl = SmilesMolSupplier('in.smi')
           >>> nMols = len(suppl)
           >>> for i in range(nMols):
           ...   suppl[i].GetNumAtoms()

      If the input file has a title line and more than two columns (smiles and id), the
      additional columns will be used to set properties on each molecule.  The properties
      are accessible using the mol.GetProp(propName) method.
    """
    def GetItemText(self, index: int) -> str: 
        """
        GetItemText( self: SmilesMolSupplier, index: int) -> str
            returns the text for an item

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::SmilesMolSupplier {lvalue},unsigned int)
        """
    def SetData(self, data: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None: 
        """
        SetData( self: SmilesMolSupplier, data: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None
            Sets the text to be parsed

            C++ signature :
                void SetData(RDKit::SmilesMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
        """
    @staticmethod
    def __enter__( arg1: SmilesMolSupplier) -> SmilesMolSupplier: 
        """
        __enter__( arg1: SmilesMolSupplier) -> SmilesMolSupplier

            C++ signature :
                RDKit::SmilesMolSupplier* __enter__(RDKit::SmilesMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: SmilesMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: SmilesMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::SmilesMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def __getitem__( arg1: SmilesMolSupplier, arg2: int) -> Mol: 
        """
        __getitem__( arg1: SmilesMolSupplier, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__(RDKit::SmilesMolSupplier*,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, data: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None: 
        """
        __init__( arg1: object, data: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> None
            Constructor
            
              ARGUMENTS: 
            
                - fileName: name of the file to be read
            
                - delimiter: (optional) text delimiter (a string).  Defauts to ' '.
            
                - smilesColumn: (optional) index of the column containing the SMILES
                  data.  Defaults to 0.
            
                - nameColumn: (optional) index of the column containing molecule names.
                  Defaults to 1.
            
                - titleLine: (optional) set this toggle if the file contains a title line.
                  Defaults to 1.
            
                - sanitize: (optional) toggles sanitization of molecules as they are read.
                  Defaults to 1.
            
            

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,int=0 [,int=1 [,bool=True [,bool=True]]]]])

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: ...
    @staticmethod
    def __iter__( arg1: SmilesMolSupplier) -> SmilesMolSupplier: 
        """
        __iter__( arg1: SmilesMolSupplier) -> SmilesMolSupplier

            C++ signature :
                RDKit::SmilesMolSupplier* __iter__(RDKit::SmilesMolSupplier*)
        """
    @staticmethod
    def __len__( arg1: SmilesMolSupplier) -> int: 
        """
        __len__( arg1: SmilesMolSupplier) -> int

            C++ signature :
                unsigned int __len__(RDKit::SmilesMolSupplier {lvalue})
        """
    @staticmethod
    def __next__( arg1: SmilesMolSupplier) -> Mol: 
        """
        __next__( arg1: SmilesMolSupplier) -> Mol
            Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::SmilesMolSupplier*)
        """
    @staticmethod
    def reset( arg1: SmilesMolSupplier) -> None: 
        """
        reset( arg1: SmilesMolSupplier) -> None
            Resets our position in the file to the beginning.
            

            C++ signature :
                void reset(RDKit::SmilesMolSupplier {lvalue})
        """
    __instance_size__ = 200
    pass
class SmilesParserParams(Boost.Python.instance):
    """
    Parameters controlling SMILES Parsing
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def allowCXSMILES(self) -> None:
        """
        controls whether or not the CXSMILES extensions are parsed

        :type: None
        """
    @property
    def debugParse(self) -> None:
        """
        controls the amount of debugging information produced

        :type: None
        """
    @property
    def parseName(self) -> None:
        """
        controls whether or not the molecule name is also parsed

        :type: None
        """
    @property
    def removeHs(self) -> None:
        """
        controls whether or not Hs are removed before the molecule is returned

        :type: None
        """
    @property
    def sanitize(self) -> None:
        """
        controls whether or not the molecule is sanitized before being returned

        :type: None
        """
    @property
    def strictCXSMILES(self) -> None:
        """
        controls whether or not problems in CXSMILES parsing causes molecule parsing to fail

        :type: None
        """
    __instance_size__ = 48
    pass
class SmilesWriteParams(Boost.Python.instance):
    """
    Parameters controlling SMILES writing
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def allBondsExplicit(self) -> None:
        """
        include symbols for all bonds

        :type: None
        """
    @property
    def allHsExplicit(self) -> None:
        """
        provide hydrogen counts for every atom

        :type: None
        """
    @property
    def canonical(self) -> None:
        """
        generate canonical SMILES

        :type: None
        """
    @property
    def doIsomericSmiles(self) -> None:
        """
        include stereochemistry and isotope information

        :type: None
        """
    @property
    def doKekule(self) -> None:
        """
        kekulize the molecule before generating the SMILES and output single/double bonds. NOTE that the output is not canonical and that this will thrown an exception if the molecule cannot be kekulized

        :type: None
        """
    @property
    def doRandom(self) -> None:
        """
        randomize the output order. The resulting SMILES is not canonical

        :type: None
        """
    @property
    def rootedAtAtom(self) -> None:
        """
        make sure the SMILES starts at the specified atom. The resulting SMILES is not canonical

        :type: None
        """
    __instance_size__ = 40
    pass
class SmilesWriter(Boost.Python.instance):
    """
    A class for writing molecules to text files.
    """
    @staticmethod
    def NumMols( arg1: SmilesWriter) -> int: 
        """
        NumMols( arg1: SmilesWriter) -> int
            Returns the number of molecules written so far.
            
            

            C++ signature :
                unsigned int NumMols(RDKit::SmilesWriter {lvalue})
        """
    @staticmethod
    def SetProps( arg1: SmilesWriter, arg2: object) -> None: 
        """
        SetProps( arg1: SmilesWriter, arg2: object) -> None
            Sets the properties to be written to the output file
            
              ARGUMENTS:
            
                - props: a list or tuple of property names
            
            

            C++ signature :
                void SetProps(RDKit::SmilesWriter {lvalue},boost::python::api::object)
        """
    @staticmethod
    def __enter__( arg1: SmilesWriter) -> SmilesWriter: 
        """
        __enter__( arg1: SmilesWriter) -> SmilesWriter

            C++ signature :
                RDKit::SmilesWriter* __enter__(RDKit::SmilesWriter*)
        """
    @staticmethod
    def __exit__( arg1: SmilesWriter, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: SmilesWriter, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::SmilesWriter*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileObj: object, delimiter: str = '', nameHeader: str = 'Name', includeHeader: bool = True, isomericSmiles: bool = True, kekuleSmiles: bool = False) -> object: 
        """
        __init__( arg1: object, fileObj: object, delimiter: str = '', nameHeader: str = 'Name', includeHeader: bool = True, isomericSmiles: bool = True, kekuleSmiles: bool = False) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object {lvalue} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='Name' [,bool=True [,bool=True [,bool=False]]]]])

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='Name' [,bool=True [,bool=True [,bool=False]]]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, delimiter: str = '', nameHeader: str = 'Name', includeHeader: bool = True, isomericSmiles: bool = True, kekuleSmiles: bool = False) -> None: ...
    @staticmethod
    def close( arg1: SmilesWriter) -> None: 
        """
        close( arg1: SmilesWriter) -> None
            Flushes the output file and closes it. The Writer cannot be used after this.
            
            

            C++ signature :
                void close(RDKit::SmilesWriter {lvalue})
        """
    @staticmethod
    def flush( arg1: SmilesWriter) -> None: 
        """
        flush( arg1: SmilesWriter) -> None
            Flushes the output file (forces the disk file to be updated).
            
            

            C++ signature :
                void flush(RDKit::SmilesWriter {lvalue})
        """
    def write(self, mol: Mol, confId: int = -1) -> None: 
        """
        write( self: SmilesWriter, mol: Mol, confId: int = -1) -> None
            Writes a molecule to the output file.
            
              ARGUMENTS:
            
                - mol: the Mol to be written
                - confId: (optional) ignored 
            
            

            C++ signature :
                void write(RDKit::SmilesWriter {lvalue},RDKit::ROMol [,int=-1])
        """
    pass
class TDTMolSupplier(Boost.Python.instance):
    """
    A class which supplies molecules from a TDT file.

      Usage examples:

        1) Lazy evaluation: the molecules are not constructed until we ask for them:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> for mol in suppl:
           ...    mol.GetNumAtoms()

        2) Lazy evaluation 2:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> mol1 = next(suppl)
           >>> mol2 = next(suppl)
           >>> suppl.reset()
           >>> mol3 = next(suppl)

           # mol3 and mol1 are the same:       >>> MolToSmiles(mol3)==MolToSmiles(mol1)

        3) Random Access:  all molecules are constructed as soon as we ask for the
           length:

           >>> suppl = TDTMolSupplier('in.smi')
           >>> nMols = len(suppl)
           >>> for i in range(nMols):
           ...   suppl[i].GetNumAtoms()

      Properties in the file are used to set properties on each molecule.
      The properties are accessible using the mol.GetProp(propName) method.
    """
    def GetItemText(self, index: int) -> str: 
        """
        GetItemText( self: TDTMolSupplier, index: int) -> str
            returns the text for an item

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetItemText(RDKit::TDTMolSupplier {lvalue},unsigned int)
        """
    def SetData(self, data: str, nameRecord: str = '', confId2D: int = -1, confId3D: int = -1, sanitize: bool = True) -> None: 
        """
        SetData( self: TDTMolSupplier, data: str, nameRecord: str = '', confId2D: int = -1, confId3D: int = -1, sanitize: bool = True) -> None
            Sets the text to be parsed

            C++ signature :
                void SetData(RDKit::TDTMolSupplier {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,int=-1 [,int=-1 [,bool=True]]]])
        """
    @staticmethod
    def __enter__( arg1: TDTMolSupplier) -> TDTMolSupplier: 
        """
        __enter__( arg1: TDTMolSupplier) -> TDTMolSupplier

            C++ signature :
                RDKit::TDTMolSupplier* __enter__(RDKit::TDTMolSupplier*)
        """
    @staticmethod
    def __exit__( arg1: TDTMolSupplier, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: TDTMolSupplier, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::TDTMolSupplier*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    def __getitem__( arg1: TDTMolSupplier, arg2: int) -> Mol: 
        """
        __getitem__( arg1: TDTMolSupplier, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__(RDKit::TDTMolSupplier*,int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,int=-1 [,int=-1 [,bool=True]]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str, nameRecord: str = '', confId2D: int = -1, confId3D: int = -1, sanitize: bool = True) -> None: ...
    @staticmethod
    def __iter__( arg1: TDTMolSupplier) -> TDTMolSupplier: 
        """
        __iter__( arg1: TDTMolSupplier) -> TDTMolSupplier

            C++ signature :
                RDKit::TDTMolSupplier* __iter__(RDKit::TDTMolSupplier*)
        """
    @staticmethod
    def __len__( arg1: TDTMolSupplier) -> int: 
        """
        __len__( arg1: TDTMolSupplier) -> int

            C++ signature :
                unsigned int __len__(RDKit::TDTMolSupplier {lvalue})
        """
    @staticmethod
    def __next__( arg1: TDTMolSupplier) -> Mol: 
        """
        __next__( arg1: TDTMolSupplier) -> Mol
            Returns the next molecule in the file.  Raises _StopIteration_ on EOF.
            

            C++ signature :
                RDKit::ROMol* __next__(RDKit::TDTMolSupplier*)
        """
    @staticmethod
    def reset( arg1: TDTMolSupplier) -> None: 
        """
        reset( arg1: TDTMolSupplier) -> None
            Resets our position in the file to the beginning.
            

            C++ signature :
                void reset(RDKit::TDTMolSupplier {lvalue})
        """
    __instance_size__ = 128
    pass
class TDTWriter(Boost.Python.instance):
    """
    A class for writing molecules to TDT files.
    """
    @staticmethod
    def GetNumDigits( arg1: TDTWriter) -> int: 
        """
        GetNumDigits( arg1: TDTWriter) -> int

            C++ signature :
                unsigned int GetNumDigits(RDKit::TDTWriter {lvalue})
        """
    @staticmethod
    def GetWrite2D( arg1: TDTWriter) -> bool: 
        """
        GetWrite2D( arg1: TDTWriter) -> bool

            C++ signature :
                bool GetWrite2D(RDKit::TDTWriter {lvalue})
        """
    @staticmethod
    def GetWriteNames( arg1: TDTWriter) -> bool: 
        """
        GetWriteNames( arg1: TDTWriter) -> bool

            C++ signature :
                bool GetWriteNames(RDKit::TDTWriter {lvalue})
        """
    @staticmethod
    def NumMols( arg1: TDTWriter) -> int: 
        """
        NumMols( arg1: TDTWriter) -> int
            Returns the number of molecules written so far.
            
            

            C++ signature :
                unsigned int NumMols(RDKit::TDTWriter {lvalue})
        """
    @staticmethod
    def SetNumDigits( arg1: TDTWriter, arg2: int) -> None: 
        """
        SetNumDigits( arg1: TDTWriter, arg2: int) -> None
            sets the number of digits to be written for coordinates

            C++ signature :
                void SetNumDigits(RDKit::TDTWriter {lvalue},unsigned int)
        """
    @staticmethod
    def SetProps( arg1: TDTWriter, arg2: object) -> None: 
        """
        SetProps( arg1: TDTWriter, arg2: object) -> None
            Sets the properties to be written to the output file
            
              ARGUMENTS:
            
                - props: a list or tuple of property names
            
            

            C++ signature :
                void SetProps(RDKit::TDTWriter {lvalue},boost::python::api::object)
        """
    def SetWrite2D(self, state: bool = True) -> None: 
        """
        SetWrite2D( self: TDTWriter, state: bool = True) -> None
            causes 2D conformations to be written (default is 3D conformations)

            C++ signature :
                void SetWrite2D(RDKit::TDTWriter {lvalue} [,bool=True])
        """
    def SetWriteNames(self, state: bool = True) -> None: 
        """
        SetWriteNames( self: TDTWriter, state: bool = True) -> None
            causes names to be written to the output file as NAME records

            C++ signature :
                void SetWriteNames(RDKit::TDTWriter {lvalue} [,bool=True])
        """
    @staticmethod
    def __enter__( arg1: TDTWriter) -> TDTWriter: 
        """
        __enter__( arg1: TDTWriter) -> TDTWriter

            C++ signature :
                RDKit::TDTWriter* __enter__(RDKit::TDTWriter*)
        """
    @staticmethod
    def __exit__( arg1: TDTWriter, arg2: object, arg3: object, arg4: object) -> bool: 
        """
        __exit__( arg1: TDTWriter, arg2: object, arg3: object, arg4: object) -> bool

            C++ signature :
                bool __exit__(RDKit::TDTWriter*,boost::python::api::object,boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: object) -> object: 
        """
        __init__( arg1: object, arg2: object) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object {lvalue})

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fileName: str) -> None: ...
    @staticmethod
    def close( arg1: TDTWriter) -> None: 
        """
        close( arg1: TDTWriter) -> None
            Flushes the output file and closes it. The Writer cannot be used after this.
            
            

            C++ signature :
                void close(RDKit::TDTWriter {lvalue})
        """
    @staticmethod
    def flush( arg1: TDTWriter) -> None: 
        """
        flush( arg1: TDTWriter) -> None
            Flushes the output file (forces the disk file to be updated).
            
            

            C++ signature :
                void flush(RDKit::TDTWriter {lvalue})
        """
    def write(self, mol: Mol, confId: int = -1) -> None: 
        """
        write( self: TDTWriter, mol: Mol, confId: int = -1) -> None
            Writes a molecule to the output file.
            
              ARGUMENTS:
            
                - mol: the Mol to be written
                - confId: (optional) ID of the conformation to write
            
            

            C++ signature :
                void write(RDKit::TDTWriter {lvalue},RDKit::ROMol [,int=-1])
        """
    pass
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
def MolToCXSmiles( mol: Mol, params: SmilesWriteParams, flags: int = CXSmilesFields.CX_ALL) -> str:
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
def SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier:
    """
    SmilesMolSupplierFromText( text: str, delimiter: str = '', smilesColumn: int = 0, nameColumn: int = 1, titleLine: bool = True, sanitize: bool = True) -> SmilesMolSupplier

        C++ signature :
            RDKit::SmilesMolSupplier* SmilesMolSupplierFromText(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=' ' [,int=0 [,int=1 [,bool=True [,bool=True]]]]])
    """
