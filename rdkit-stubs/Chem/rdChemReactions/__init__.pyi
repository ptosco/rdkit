"""Module containing classes and functions for working with chemical reactions."""
from __future__ import annotations
import rdkit.Chem.rdChemReactions
import typing
import Boost.Python

__all__ = [
    "CartesianProductStrategy",
    "ChemicalReaction",
    "Compute2DCoordsForReaction",
    "CreateDifferenceFingerprintForReaction",
    "CreateStructuralFingerprintForReaction",
    "EnumerateLibrary",
    "EnumerateLibraryBase",
    "EnumerateLibraryCanSerialize",
    "EnumerationParams",
    "EnumerationStrategyBase",
    "EvenSamplePairsStrategy",
    "FingerprintType",
    "GetChemDrawRxnAdjustParams",
    "GetDefaultAdjustParams",
    "HasAgentTemplateSubstructMatch",
    "HasProductTemplateSubstructMatch",
    "HasReactantTemplateSubstructMatch",
    "HasReactionAtomMapping",
    "HasReactionSubstructMatch",
    "IsReactionTemplateMoleculeAgent",
    "MOL_SPTR_VECT",
    "MatchOnlyAtRgroupsAdjustParams",
    "MrvBlockIsReaction",
    "MrvFileIsReaction",
    "PreprocessReaction",
    "RandomSampleAllBBsStrategy",
    "RandomSampleStrategy",
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
    "ReduceProductToSideChains",
    "RemoveMappingNumbersFromReactions",
    "SANITIZE_ADJUST_REACTANTS",
    "SANITIZE_ALL",
    "SANITIZE_ATOM_MAPS",
    "SANITIZE_MERGEHS",
    "SANITIZE_NONE",
    "SANITIZE_RGROUP_NAMES",
    "SanitizeFlags",
    "SanitizeRxn",
    "UpdateProductsStereochemistry",
    "VectMolVect",
    "VectSizeT",
    "VectorOfStringVectors"
]


class EnumerationStrategyBase(Boost.Python.instance):
    @staticmethod
    def GetNumPermutations( arg1: EnumerationStrategyBase) -> int: 
        """
        GetNumPermutations( arg1: EnumerationStrategyBase) -> int
            Returns the total number of results for this enumeration strategy.
            Note that some strategies are effectively infinite.

            C++ signature :
                unsigned long GetNumPermutations(RDKit::EnumerationStrategyBase {lvalue})
        """
    @staticmethod
    def GetPosition( arg1: EnumerationStrategyBase) -> VectSizeT: 
        """
        GetPosition( arg1: EnumerationStrategyBase) -> VectSizeT
            Return the current indices into the arrays of reagents

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > GetPosition(RDKit::EnumerationStrategyBase {lvalue})
        """
    @staticmethod
    def Initialize( arg1: EnumerationStrategyBase, arg2: ChemicalReaction, arg3: list) -> None: 
        """
        Initialize( arg1: EnumerationStrategyBase, arg2: ChemicalReaction, arg3: list) -> None

            C++ signature :
                void Initialize(RDKit::EnumerationStrategyBase {lvalue},RDKit::ChemicalReaction {lvalue},boost::python::list)
        """
    @staticmethod
    def Skip( arg1: EnumerationStrategyBase, skipCount: int) -> bool: 
        """
        Skip( arg1: EnumerationStrategyBase, skipCount: int) -> bool
            Skip the next Nth results. note: this may be an expensive operation
            depending on the enumeration strategy used. It is recommended to use
            the enumerator state to advance to a known position

            C++ signature :
                bool Skip(RDKit::EnumerationStrategyBase {lvalue},unsigned long)
        """
    @staticmethod
    def Type( arg1: EnumerationStrategyBase) -> str: 
        """
        Type( arg1: EnumerationStrategyBase) -> str
            Returns the enumeration strategy type as a string.

            C++ signature :
                char const* Type(RDKit::EnumerationStrategyBase {lvalue})
        """
    @staticmethod
    def __bool__( arg1: EnumerationStrategyBase) -> bool: 
        """
        __bool__( arg1: EnumerationStrategyBase) -> bool

            C++ signature :
                bool __bool__(RDKit::EnumerationStrategyBase*)
        """
    @staticmethod
    @typing.overload
    def __copy__( arg1: EnumerationStrategyBase) -> EnumerationStrategyBase: 
        """
        __copy__( arg1: EnumerationStrategyBase) -> EnumerationStrategyBase

            C++ signature :
                RDKit::EnumerationStrategyBase* __copy__(RDKit::EnumerationStrategyBase {lvalue})

            C++ signature :
                void __copy__(RDKit::EnumerationStrategyBase* {lvalue})
        """
    @staticmethod
    @typing.overload
    def __copy__( arg1: EnumerationStrategyBase) -> None: ...
    @staticmethod
    @typing.overload
    def __next__( arg1: EnumerationStrategyBase) -> VectSizeT: 
        """
        __next__( arg1: EnumerationStrategyBase) -> VectSizeT
            Return the next indices into the arrays of reagents

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > __next__(RDKit::EnumerationStrategyBase {lvalue})

            C++ signature :
                void __next__(RDKit::EnumerationStrategyBase* {lvalue})
        """
    @staticmethod
    @typing.overload
    def __next__( arg1: EnumerationStrategyBase) -> None: ...
    @staticmethod
    def __nonzero__( arg1: EnumerationStrategyBase) -> bool: 
        """
        __nonzero__( arg1: EnumerationStrategyBase) -> bool

            C++ signature :
                bool __nonzero__(RDKit::EnumerationStrategyBase*)
        """
    @staticmethod
    @typing.overload
    def next( arg1: EnumerationStrategyBase) -> VectSizeT: 
        """
        next( arg1: EnumerationStrategyBase) -> VectSizeT
            Return the next indices into the arrays of reagents

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > next(RDKit::EnumerationStrategyBase {lvalue})

            C++ signature :
                void next(RDKit::EnumerationStrategyBase* {lvalue})
        """
    @staticmethod
    @typing.overload
    def next( arg1: EnumerationStrategyBase) -> None: ...
    pass
class ChemicalReaction(Boost.Python.instance):
    """
    A class for storing and applying chemical reactions.

    Sample Usage:
      >>> from rdkit import Chem
      >>> from rdkit.Chem import rdChemReactions
      >>> rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')
      >>> reacts = (Chem.MolFromSmiles('C(=O)O'),Chem.MolFromSmiles('CNC'))
      >>> products = rxn.RunReactants(reacts)
      >>> len(products)
      1
      >>> len(products[0])
      1
      >>> Chem.MolToSmiles(products[0][0])
      'CN(C)C=O'
    """
    @staticmethod
    def AddAgentTemplate( arg1: ChemicalReaction, arg2: Mol) -> int: 
        """
        AddAgentTemplate( arg1: ChemicalReaction, arg2: Mol) -> int
            adds a agent (a Molecule)

            C++ signature :
                unsigned int AddAgentTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def AddProductTemplate( arg1: ChemicalReaction, arg2: Mol) -> int: 
        """
        AddProductTemplate( arg1: ChemicalReaction, arg2: Mol) -> int
            adds a product (a Molecule)

            C++ signature :
                unsigned int AddProductTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def AddReactantTemplate( arg1: ChemicalReaction, arg2: Mol) -> int: 
        """
        AddReactantTemplate( arg1: ChemicalReaction, arg2: Mol) -> int
            adds a reactant (a Molecule) to the reaction

            C++ signature :
                unsigned int AddReactantTemplate(RDKit::ChemicalReaction {lvalue},boost::shared_ptr<RDKit::ROMol>)
        """
    @staticmethod
    def AddRecursiveQueriesToReaction( reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue', getLabels: bool = False) -> object: 
        """
        AddRecursiveQueriesToReaction( reaction: ChemicalReaction, queries: dict = {}, propName: str = 'molFileValue', getLabels: bool = False) -> object
            adds recursive queries and returns reactant labels

            C++ signature :
                boost::python::api::object AddRecursiveQueriesToReaction(RDKit::ChemicalReaction {lvalue} [,boost::python::dict={} [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='molFileValue' [,bool=False]]])
        """
    @staticmethod
    def ClearComputedProps( arg1: ChemicalReaction) -> None: 
        """
        ClearComputedProps( arg1: ChemicalReaction) -> None
            Removes all computed properties from the reaction.
            
            

            C++ signature :
                void ClearComputedProps(RDKit::ChemicalReaction)
        """
    @staticmethod
    def ClearProp( arg1: ChemicalReaction, arg2: str) -> None: 
        """
        ClearProp( arg1: ChemicalReaction, arg2: str) -> None
            Removes a property from the reaction.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            

            C++ signature :
                void ClearProp(RDKit::ChemicalReaction,char const*)
        """
    def GetAgentTemplate(self, which: int) -> Mol: 
        """
        GetAgentTemplate( self: ChemicalReaction, which: int) -> Mol
            returns one of our agent templates

            C++ signature :
                RDKit::ROMol* GetAgentTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
    @staticmethod
    def GetAgents( arg1: ChemicalReaction) -> MOL_SPTR_VECT: 
        """
        GetAgents( arg1: ChemicalReaction) -> MOL_SPTR_VECT
            get the agent templates

            C++ signature :
                std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetAgents(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def GetBoolProp( arg1: ChemicalReaction, arg2: str) -> bool: 
        """
        GetBoolProp( arg1: ChemicalReaction, arg2: str) -> bool
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                bool GetBoolProp(RDKit::ChemicalReaction const*,char const*)
        """
    @staticmethod
    def GetDoubleProp( arg1: ChemicalReaction, arg2: str) -> float: 
        """
        GetDoubleProp( arg1: ChemicalReaction, arg2: str) -> float
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                double GetDoubleProp(RDKit::ChemicalReaction const*,char const*)
        """
    @staticmethod
    def GetIntProp( arg1: ChemicalReaction, arg2: str) -> int: 
        """
        GetIntProp( arg1: ChemicalReaction, arg2: str) -> int
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                int GetIntProp(RDKit::ChemicalReaction const*,char const*)
        """
    @staticmethod
    def GetNumAgentTemplates( arg1: ChemicalReaction) -> int: 
        """
        GetNumAgentTemplates( arg1: ChemicalReaction) -> int
            returns the number of agents this reaction expects

            C++ signature :
                unsigned int GetNumAgentTemplates(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def GetNumProductTemplates( arg1: ChemicalReaction) -> int: 
        """
        GetNumProductTemplates( arg1: ChemicalReaction) -> int
            returns the number of products this reaction generates

            C++ signature :
                unsigned int GetNumProductTemplates(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def GetNumReactantTemplates( arg1: ChemicalReaction) -> int: 
        """
        GetNumReactantTemplates( arg1: ChemicalReaction) -> int
            returns the number of reactants this reaction expects

            C++ signature :
                unsigned int GetNumReactantTemplates(RDKit::ChemicalReaction {lvalue})
        """
    def GetProductTemplate(self, which: int) -> Mol: 
        """
        GetProductTemplate( self: ChemicalReaction, which: int) -> Mol
            returns one of our product templates

            C++ signature :
                RDKit::ROMol* GetProductTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
    @staticmethod
    def GetProducts( arg1: ChemicalReaction) -> MOL_SPTR_VECT: 
        """
        GetProducts( arg1: ChemicalReaction) -> MOL_SPTR_VECT
            get the product templates

            C++ signature :
                std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetProducts(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def GetProp( arg1: ChemicalReaction, arg2: str) -> str: 
        """
        GetProp( arg1: ChemicalReaction, arg2: str) -> str
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetProp(RDKit::ChemicalReaction const*,char const*)
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE: 
        """
        GetPropNames( self: ChemicalReaction, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE
            Returns a tuple with all property names for this reaction.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            

            C++ signature :
                std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > GetPropNames(RDKit::ChemicalReaction {lvalue} [,bool=False [,bool=False]])
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict: 
        """
        GetPropsAsDict( self: ChemicalReaction, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict
            Returns a dictionary populated with the reaction's properties.
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
              RETURNS: a dictionary
            

            C++ signature :
                boost::python::dict GetPropsAsDict(RDKit::ChemicalReaction [,bool=False [,bool=False [,bool=True]]])
        """
    def GetReactantTemplate(self, which: int) -> Mol: 
        """
        GetReactantTemplate( self: ChemicalReaction, which: int) -> Mol
            returns one of our reactant templates

            C++ signature :
                RDKit::ROMol* GetReactantTemplate(RDKit::ChemicalReaction const*,unsigned int)
        """
    @staticmethod
    def GetReactants( arg1: ChemicalReaction) -> MOL_SPTR_VECT: 
        """
        GetReactants( arg1: ChemicalReaction) -> MOL_SPTR_VECT
            get the reactant templates

            C++ signature :
                std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > GetReactants(RDKit::ChemicalReaction {lvalue})
        """
    def GetReactingAtoms(self, mappedAtomsOnly: bool = False) -> object: 
        """
        GetReactingAtoms( self: ChemicalReaction, mappedAtomsOnly: bool = False) -> object
            returns a sequence of sequences with the atoms that change in the reaction

            C++ signature :
                boost::python::api::object GetReactingAtoms(RDKit::ChemicalReaction [,bool=False])
        """
    @staticmethod
    def GetSubstructParams( arg1: ChemicalReaction) -> SubstructMatchParameters: 
        """
        GetSubstructParams( arg1: ChemicalReaction) -> SubstructMatchParameters
            get the parameter object controlling the substructure matching

            C++ signature :
                RDKit::SubstructMatchParameters* GetSubstructParams(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def GetUnsignedProp( arg1: ChemicalReaction, arg2: str) -> int: 
        """
        GetUnsignedProp( arg1: ChemicalReaction, arg2: str) -> int
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            

            C++ signature :
                unsigned int GetUnsignedProp(RDKit::ChemicalReaction const*,char const*)
        """
    @staticmethod
    def HasProp( arg1: ChemicalReaction, arg2: str) -> int: 
        """
        HasProp( arg1: ChemicalReaction, arg2: str) -> int
            Queries a molecule to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            

            C++ signature :
                int HasProp(RDKit::ChemicalReaction,char const*)
        """
    def Initialize(self, silent: bool = False) -> None: 
        """
        Initialize( self: ChemicalReaction, silent: bool = False) -> None
            initializes the reaction so that it can be used

            C++ signature :
                void Initialize(RDKit::ChemicalReaction {lvalue} [,bool=False])
        """
    @staticmethod
    def IsInitialized( arg1: ChemicalReaction) -> bool: 
        """
        IsInitialized( arg1: ChemicalReaction) -> bool
            checks if the reaction is ready for use

            C++ signature :
                bool IsInitialized(RDKit::ChemicalReaction {lvalue})
        """
    @staticmethod
    def IsMoleculeAgent( arg1: ChemicalReaction, arg2: Mol) -> bool: 
        """
        IsMoleculeAgent( arg1: ChemicalReaction, arg2: Mol) -> bool
            returns whether or not the molecule has a substructure match to one of the agents.

            C++ signature :
                bool IsMoleculeAgent(RDKit::ChemicalReaction,RDKit::ROMol)
        """
    @staticmethod
    def IsMoleculeProduct( arg1: ChemicalReaction, arg2: Mol) -> bool: 
        """
        IsMoleculeProduct( arg1: ChemicalReaction, arg2: Mol) -> bool
            returns whether or not the molecule has a substructure match to one of the products.

            C++ signature :
                bool IsMoleculeProduct(RDKit::ChemicalReaction,RDKit::ROMol)
        """
    @staticmethod
    def IsMoleculeReactant( arg1: ChemicalReaction, arg2: Mol) -> bool: 
        """
        IsMoleculeReactant( arg1: ChemicalReaction, arg2: Mol) -> bool
            returns whether or not the molecule has a substructure match to one of the reactants.

            C++ signature :
                bool IsMoleculeReactant(RDKit::ChemicalReaction,RDKit::ROMol)
        """
    def RemoveAgentTemplates(self, targetList: AtomPairsParameters = None) -> None: 
        """
        RemoveAgentTemplates( self: ChemicalReaction, targetList: AtomPairsParameters = None) -> None
            Removes agents from reaction. If targetList is provide the agents will be transferred to that list.

            C++ signature :
                void RemoveAgentTemplates(RDKit::ChemicalReaction {lvalue} [,boost::python::api::object=None])
        """
    def RemoveUnmappedProductTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: AtomPairsParameters = None) -> None: 
        """
        RemoveUnmappedProductTemplates( self: ChemicalReaction, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: AtomPairsParameters = None) -> None
            Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from product templates to the agent templates or to a given targetList

            C++ signature :
                void RemoveUnmappedProductTemplates(RDKit::ChemicalReaction* [,double=0.2 [,bool=True [,boost::python::api::object=None]]])
        """
    def RemoveUnmappedReactantTemplates(self, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: AtomPairsParameters = None) -> None: 
        """
        RemoveUnmappedReactantTemplates( self: ChemicalReaction, thresholdUnmappedAtoms: float = 0.2, moveToAgentTemplates: bool = True, targetList: AtomPairsParameters = None) -> None
            Removes molecules with an atom mapping ratio below thresholdUnmappedAtoms from reactant templates to the agent templates or to a given targetList

            C++ signature :
                void RemoveUnmappedReactantTemplates(RDKit::ChemicalReaction* [,double=0.2 [,bool=True [,boost::python::api::object=None]]])
        """
    @staticmethod
    def RunReactant( arg1: ChemicalReaction, arg2: AtomPairsParameters, arg3: int) -> object: 
        """
        RunReactant( arg1: ChemicalReaction, arg2: AtomPairsParameters, arg3: int) -> object
            apply the reaction to a single reactant

            C++ signature :
                _object* RunReactant(RDKit::ChemicalReaction*,boost::python::api::object,unsigned int)
        """
    def RunReactantInPlace(self, reactant: Mol, removeUnmatchedAtoms: bool = True) -> bool: 
        """
        RunReactantInPlace( self: ChemicalReaction, reactant: Mol, removeUnmatchedAtoms: bool = True) -> bool
            apply the reaction to a single reactant in place. The reactant itself is modified. This can only be used for single reactant - single product reactions.

            C++ signature :
                bool RunReactantInPlace(RDKit::ChemicalReaction*,RDKit::ROMol* [,bool=True])
        """
    @typing.overload
    def RunReactants(self, reactants: tuple, maxProducts: int = 1000) -> object: 
        """
        RunReactants( self: ChemicalReaction, reactants: tuple, maxProducts: int = 1000) -> object
            apply the reaction to a sequence of reactant molecules and return the products as a tuple of tuples.  If maxProducts is not zero, stop the reaction when maxProducts have been generated [default=1000]

            C++ signature :
                _object* RunReactants(RDKit::ChemicalReaction*,boost::python::tuple [,unsigned int=1000])

            C++ signature :
                _object* RunReactants(RDKit::ChemicalReaction*,boost::python::list [,unsigned int=1000])
        """
    @typing.overload
    def RunReactants(self, reactants: list, maxProducts: int = 1000) -> object: ...
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None: 
        """
        SetBoolProp( self: ChemicalReaction, key: str, val: bool, computed: bool = False) -> None
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetBoolProp(RDKit::ChemicalReaction,char const*,bool [,bool=False])
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None: 
        """
        SetDoubleProp( self: ChemicalReaction, key: str, val: float, computed: bool = False) -> None
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            

            C++ signature :
                void SetDoubleProp(RDKit::ChemicalReaction,char const*,double [,bool=False])
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetIntProp( self: ChemicalReaction, key: str, val: int, computed: bool = False) -> None
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetIntProp(RDKit::ChemicalReaction,char const*,int [,bool=False])
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None: 
        """
        SetProp( self: ChemicalReaction, key: str, val: str, computed: bool = False) -> None
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetProp(RDKit::ChemicalReaction,char const*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None: 
        """
        SetUnsignedProp( self: ChemicalReaction, key: str, val: int, computed: bool = False) -> None
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            

            C++ signature :
                void SetUnsignedProp(RDKit::ChemicalReaction,char const*,unsigned int [,bool=False])
        """
    @typing.overload
    def ToBinary(self) -> object: 
        """
        ToBinary( self: ChemicalReaction) -> object
            Returns a binary string representation of the reaction.

            C++ signature :
                boost::python::api::object ToBinary(RDKit::ChemicalReaction)

            C++ signature :
                boost::python::api::object ToBinary(RDKit::ChemicalReaction,unsigned int)
        """
    @typing.overload
    def ToBinary(self, propertyFlags: int) -> object: ...
    def Validate(self, silent: bool = False) -> tuple: 
        """
        Validate( self: ChemicalReaction, silent: bool = False) -> tuple
            checks the reaction for potential problems, returns (numWarnings,numErrors)

            C++ signature :
                boost::python::tuple Validate(RDKit::ChemicalReaction const* [,bool=False])
        """
    @staticmethod
    def __getinitargs__( arg1: ChemicalReaction) -> tuple: 
        """
        __getinitargs__( arg1: ChemicalReaction) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::ChemicalReaction)
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
            Constructor, takes no arguments

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,RDKit::ChemicalReaction)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: ChemicalReaction) -> None: ...
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    def _getImplicitPropertiesFlag(self) -> bool: 
        """
        _getImplicitPropertiesFlag( self: ChemicalReaction) -> bool
            EXPERT USER: returns whether or not the reaction can have implicit properties

            C++ signature :
                bool _getImplicitPropertiesFlag(RDKit::ChemicalReaction {lvalue})
        """
    def _setImplicitPropertiesFlag(self, val: bool) -> None: 
        """
        _setImplicitPropertiesFlag( self: ChemicalReaction, val: bool) -> None
            EXPERT USER: indicates that the reaction can have implicit properties

            C++ signature :
                void _setImplicitPropertiesFlag(RDKit::ChemicalReaction {lvalue},bool)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 40
    __safe_for_unpickling__ = True
    pass
class EnumerateLibraryBase(Boost.Python.instance):
    @staticmethod
    def GetEnumerator( arg1: EnumerateLibraryBase) -> EnumerationStrategyBase: 
        """
        GetEnumerator( arg1: EnumerateLibraryBase) -> EnumerationStrategyBase
            Returns the enumation strategy for the current library

            C++ signature :
                RDKit::EnumerationStrategyBase GetEnumerator(RDKit::EnumerateLibraryBase {lvalue})
        """
    @staticmethod
    def GetPosition( arg1: EnumerateLibraryBase) -> VectSizeT: 
        """
        GetPosition( arg1: EnumerateLibraryBase) -> VectSizeT
            Returns the current enumeration position into the reagent vectors

            C++ signature :
                std::vector<unsigned long, std::allocator<unsigned long> > GetPosition(RDKit::EnumerateLibraryBase {lvalue})
        """
    @staticmethod
    def GetReaction( arg1: EnumerateLibraryBase) -> ChemicalReaction: 
        """
        GetReaction( arg1: EnumerateLibraryBase) -> ChemicalReaction
            Returns the chemical reaction for this library

            C++ signature :
                RDKit::ChemicalReaction GetReaction(RDKit::EnumerateLibraryBase {lvalue})
        """
    @staticmethod
    def GetState( arg1: EnumerateLibraryBase) -> str: 
        """
        GetState( arg1: EnumerateLibraryBase) -> str
            Returns the current enumeration state (position) of the library.
            This position can be used to restart the library from a known position

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetState(RDKit::EnumerateLibraryBase {lvalue})
        """
    @staticmethod
    def InitFromString( arg1: EnumerateLibraryBase, data: str) -> None: 
        """
        InitFromString( arg1: EnumerateLibraryBase, data: str) -> None
            Inititialize the library from a binary string

            C++ signature :
                void InitFromString(RDKit::EnumerateLibraryBase {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def ResetState( arg1: EnumerateLibraryBase) -> None: 
        """
        ResetState( arg1: EnumerateLibraryBase) -> None
            Returns the current enumeration state (position) of the library to the start.

            C++ signature :
                void ResetState(RDKit::EnumerateLibraryBase {lvalue})
        """
    @staticmethod
    def Serialize( arg1: EnumerateLibraryBase) -> object: 
        """
        Serialize( arg1: EnumerateLibraryBase) -> object
            Serialize the library to a binary string.
            Note that the position in the library is serialized as well.  Care should
            be taken when serializing.  See GetState/SetState for position manipulation.

            C++ signature :
                boost::python::api::object Serialize(RDKit::EnumerateLibraryBase)
        """
    @staticmethod
    def SetState( arg1: EnumerateLibraryBase, state: str) -> None: 
        """
        SetState( arg1: EnumerateLibraryBase, state: str) -> None
            Sets the enumeration state (position) of the library.

            C++ signature :
                void SetState(RDKit::EnumerateLibraryBase {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __bool__( arg1: EnumerateLibraryBase) -> bool: 
        """
        __bool__( arg1: EnumerateLibraryBase) -> bool

            C++ signature :
                bool __bool__(RDKit::EnumerateLibraryBase*)
        """
    @staticmethod
    def __iter__( arg1: AtomPairsParameters) -> object: 
        """
        __iter__( arg1: AtomPairsParameters) -> object

            C++ signature :
                boost::python::api::object __iter__(boost::python::api::object)
        """
    @staticmethod
    def __next__( arg1: EnumerateLibraryBase) -> object: 
        """
        __next__( arg1: EnumerateLibraryBase) -> object
            Return the next molecule from the enumeration.

            C++ signature :
                _object* __next__(RDKit::EnumerateLibraryBase*)
        """
    @staticmethod
    def __nonzero__( arg1: EnumerateLibraryBase) -> bool: 
        """
        __nonzero__( arg1: EnumerateLibraryBase) -> bool

            C++ signature :
                bool __nonzero__(RDKit::EnumerateLibraryBase*)
        """
    @staticmethod
    def next( arg1: EnumerateLibraryBase) -> object: 
        """
        next( arg1: EnumerateLibraryBase) -> object
            Return the next molecule from the enumeration.

            C++ signature :
                _object* next(RDKit::EnumerateLibraryBase*)
        """
    @staticmethod
    def nextSmiles( arg1: EnumerateLibraryBase) -> VectorOfStringVectors: 
        """
        nextSmiles( arg1: EnumerateLibraryBase) -> VectorOfStringVectors
            Return the next smiles string from the enumeration.

            C++ signature :
                std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > nextSmiles(RDKit::EnumerateLibraryBase {lvalue})
        """
    pass
class EnumerateLibrary(EnumerateLibraryBase, Boost.Python.instance):
    """
    EnumerateLibrary
    This class allows easy enumeration of reactions.  Simply provide a reaction
    and a set of reagents and you are off the races.

    Note that this functionality should be considered beta and that the API may
    change in a future release.

    EnumerateLibrary follows the python enumerator protocol, for example:

    library = EnumerateLibrary(rxn, bbs)
    for products in library:
       ... do something with the product

    It is useful to sanitize reactions before hand:

    SanitizeRxn(rxn)
    library = EnumerateLibrary(rxn, bbs)

    If ChemDraw style reaction semantics are prefereed, you can apply
    the ChemDraw parameters:

    SanitizeRxn(rxn, params=GetChemDrawRxnAdjustParams())

    For one, this enforces only matching RGroups and assumes all atoms
    have fully satisfied valences.

    Each product has the same output as applying a set of reagents to
    the libraries reaction.

    This can be a bit confusing as each product can have multiple molecules
    generated.  The returned data structure is as follows:

       [ [products1], [products2],... ]
    Where products1 are the molecule products for the reactions first product
    template and products2 are the molecule products for the second product
    template.  Since each reactant can match more than once, there may be
    multiple product molecules for each template.

    for products in library:
        for results_for_product_template in products:
            for mol in results_for_product_template:
                Chem.MolToSmiles(mol) # finally have a molecule!

    For sufficiently large libraries, using this iteration strategy is not
    recommended as the library may contain more products than atoms in the
    universe.  To help with this, you can supply an enumeration strategy.
    The default strategy is a CartesianProductStrategy which enumerates
    everything.  RandomSampleStrategy randomly samples the products but
    this strategy never terminates, however, python supplies itertools:

    import itertools
    library = EnumerateLibrary(rxn, bbs, rdChemReactions.RandomSampleStrategy())
    for result in itertools.islice(library, 1000):
        # do something with the first 1000 samples

    for result in itertools.islice(library, 1000):
        # do something with the next 1000 samples

    Libraries are also serializable, including their current state:

    s = library.Serialize()
    library2 = EnumerateLibrary()
    library2.InitFromString(s)
    for result in itertools.islice(libary2, 1000):
        # do something with the next 1000 samples
    """
    @staticmethod
    def GetReagents( arg1: EnumerateLibrary) -> VectMolVect: 
        """
        GetReagents( arg1: EnumerateLibrary) -> VectMolVect
            Return the reagents used in this library.

            C++ signature :
                std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > GetReagents(RDKit::EnumerateLibraryWrap {lvalue})
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,RDKit::ChemicalReaction,boost::python::list [,RDKit::EnumerationParams])

            C++ signature :
                void __init__(_object*,RDKit::ChemicalReaction,boost::python::tuple [,RDKit::EnumerationParams])

            C++ signature :
                void __init__(_object*,RDKit::ChemicalReaction,boost::python::list,RDKit::EnumerationStrategyBase [,RDKit::EnumerationParams])

            C++ signature :
                void __init__(_object*,RDKit::ChemicalReaction,boost::python::tuple,RDKit::EnumerationStrategyBase [,RDKit::EnumerationParams])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, rxn: ChemicalReaction, reagents: list, params: EnumerationParams) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, rxn: ChemicalReaction, reagents: tuple, params: EnumerationParams) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, rxn: ChemicalReaction, reagents: list, enumerator: EnumerationStrategyBase, params: EnumerationParams) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, rxn: ChemicalReaction, reagents: tuple, enumerator: EnumerationStrategyBase, params: EnumerationParams) -> None: ...
    __instance_size__ = 296
    pass
class EnumerationParams(Boost.Python.instance):
    """
    EnumerationParams
    Controls some aspects of how the enumeration is performed.
    Options:
      reagentMaxMatchCount [ default Infinite ]
        This specifies how many times the reactant template can match a reagent.

      sanePartialProducts [default false]
        If true, forces all products of the reagent plus the product templates
         pass chemical sanitization.  Note that if the product template itself
         does not pass sanitization, then none of the products will.
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def reagentMaxMatchCount(self) -> None:
        """
        :type: None
        """
    @property
    def sanePartialProducts(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 32
    pass
class CartesianProductStrategy(EnumerationStrategyBase, Boost.Python.instance):
    """
    CartesianProductStrategy produces a standard walk through all possible
    reagent combinations:

    (0,0,0), (1,0,0), (2,0,0) ...
    """
    @staticmethod
    def __copy__( arg1: CartesianProductStrategy) -> EnumerationStrategyBase: 
        """
        __copy__( arg1: CartesianProductStrategy) -> EnumerationStrategyBase

            C++ signature :
                RDKit::EnumerationStrategyBase* __copy__(RDKit::CartesianProductStrategy {lvalue})
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
class EvenSamplePairsStrategy(EnumerationStrategyBase, Boost.Python.instance):
    """
    Randomly sample Pairs evenly from a collection of building blocks
    This is a good strategy for choosing a relatively small selection
    of building blocks from a larger set.  As the amount of work needed
    to retrieve the next evenly sample building block grows with the
    number of samples, this method performs progressively worse as the
    number of samples gets larger.
    See EnumerationStrategyBase for more details.
    """
    @staticmethod
    def Stats( arg1: EvenSamplePairsStrategy) -> str: 
        """
        Stats( arg1: EvenSamplePairsStrategy) -> str
            Return the statistics log of the pairs used in the current enumeration.

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Stats(RDKit::EvenSamplePairsStrategy {lvalue})
        """
    @staticmethod
    def __copy__( arg1: EvenSamplePairsStrategy) -> EnumerationStrategyBase: 
        """
        __copy__( arg1: EvenSamplePairsStrategy) -> EnumerationStrategyBase

            C++ signature :
                RDKit::EnumerationStrategyBase* __copy__(RDKit::EvenSamplePairsStrategy {lvalue})
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
class FingerprintType(Boost.Python.enum, int):
    AtomPairFP = rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP
    MorganFP = rdkit.Chem.rdChemReactions.FingerprintType.MorganFP
    PatternFP = rdkit.Chem.rdChemReactions.FingerprintType.PatternFP
    RDKitFP = rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP
    TopologicalTorsion = rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion
    __slots__ = ()
    names = {'AtomPairFP': rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP, 'TopologicalTorsion': rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion, 'MorganFP': rdkit.Chem.rdChemReactions.FingerprintType.MorganFP, 'RDKitFP': rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP, 'PatternFP': rdkit.Chem.rdChemReactions.FingerprintType.PatternFP}
    values = {1: rdkit.Chem.rdChemReactions.FingerprintType.AtomPairFP, 2: rdkit.Chem.rdChemReactions.FingerprintType.TopologicalTorsion, 3: rdkit.Chem.rdChemReactions.FingerprintType.MorganFP, 4: rdkit.Chem.rdChemReactions.FingerprintType.RDKitFP, 5: rdkit.Chem.rdChemReactions.FingerprintType.PatternFP}
    pass
class MOL_SPTR_VECT(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: MOL_SPTR_VECT, arg2: object) -> bool: 
        """
        __contains__( arg1: MOL_SPTR_VECT, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: MOL_SPTR_VECT, arg2: object) -> None: 
        """
        __delitem__( arg1: MOL_SPTR_VECT, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<boost::shared_ptr<RDKit::ROMol>*, std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > __iter__(boost::python::back_reference<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >&>)
        """
    @staticmethod
    def __len__( arg1: MOL_SPTR_VECT) -> int: 
        """
        __len__( arg1: MOL_SPTR_VECT) -> int

            C++ signature :
                unsigned long __len__(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: MOL_SPTR_VECT, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: MOL_SPTR_VECT, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: MOL_SPTR_VECT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class RandomSampleAllBBsStrategy(EnumerationStrategyBase, Boost.Python.instance):
    """
    RandomSampleAllBBsStrategy randomly samples from the reagent sets
    with the constraint that all building blocks are samples as early as possible.
    Note that this strategy never halts and can produce duplicates.
    """
    @staticmethod
    def __copy__( arg1: RandomSampleAllBBsStrategy) -> EnumerationStrategyBase: 
        """
        __copy__( arg1: RandomSampleAllBBsStrategy) -> EnumerationStrategyBase

            C++ signature :
                RDKit::EnumerationStrategyBase* __copy__(RDKit::RandomSampleAllBBsStrategy {lvalue})
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
class RandomSampleStrategy(EnumerationStrategyBase, Boost.Python.instance):
    """
    RandomSampleStrategy simply randomly samples from the reagent sets.
    Note that this strategy never halts and can produce duplicates.
    """
    @staticmethod
    def __copy__( arg1: RandomSampleStrategy) -> EnumerationStrategyBase: 
        """
        __copy__( arg1: RandomSampleStrategy) -> EnumerationStrategyBase

            C++ signature :
                RDKit::EnumerationStrategyBase* __copy__(RDKit::RandomSampleStrategy {lvalue})
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 32
    pass
class ReactionFingerprintParams(Boost.Python.instance):
    """
    A class for storing parameters to manipulate the calculation of fingerprints of chemical reactions.
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
                void __init__(_object*,bool,double,unsigned int,int,unsigned int,RDKit::FingerprintType)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: bool, arg3: float, arg4: int, arg5: int, arg6: int, arg7: FingerprintType) -> None: ...
    @property
    def agentWeight(self) -> None:
        """
        :type: None
        """
    @property
    def bitRatioAgents(self) -> None:
        """
        :type: None
        """
    @property
    def fpSize(self) -> None:
        """
        :type: None
        """
    @property
    def fpType(self) -> None:
        """
        :type: None
        """
    @property
    def includeAgents(self) -> None:
        """
        :type: None
        """
    @property
    def nonAgentWeight(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 56
    pass
class SanitizeFlags(Boost.Python.enum, int):
    SANITIZE_ADJUST_REACTANTS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
    SANITIZE_ALL = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
    SANITIZE_ATOM_MAPS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
    SANITIZE_MERGEHS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
    SANITIZE_NONE = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
    SANITIZE_RGROUP_NAMES = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
    __slots__ = ()
    names = {'SANITIZE_NONE': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE, 'SANITIZE_ATOM_MAPS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS, 'SANITIZE_RGROUP_NAMES': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES, 'SANITIZE_ADJUST_REACTANTS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS, 'SANITIZE_MERGEHS': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS, 'SANITIZE_ALL': rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL}
    values = {0: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE, 2: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS, 1: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES, 4: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS, 8: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS, 4294967295: rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL}
    pass
class VectMolVect(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: VectMolVect, arg2: object) -> bool: 
        """
        __contains__( arg1: VectMolVect, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: VectMolVect, arg2: object) -> None: 
        """
        __delitem__( arg1: VectMolVect, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >*, std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > > > __iter__(boost::python::back_reference<std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > >&>)
        """
    @staticmethod
    def __len__( arg1: VectMolVect) -> int: 
        """
        __len__( arg1: VectMolVect) -> int

            C++ signature :
                unsigned long __len__(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: VectMolVect, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: VectMolVect, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: VectMolVect, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: VectMolVect, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: VectMolVect, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: VectMolVect, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > > > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class VectSizeT(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: VectSizeT, arg2: object) -> bool: 
        """
        __contains__( arg1: VectSizeT, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: VectSizeT, arg2: object) -> None: 
        """
        __delitem__( arg1: VectSizeT, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<unsigned long, std::allocator<unsigned long> >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<unsigned long*, std::vector<unsigned long, std::allocator<unsigned long> > > > __iter__(boost::python::back_reference<std::vector<unsigned long, std::allocator<unsigned long> >&>)
        """
    @staticmethod
    def __len__( arg1: VectSizeT) -> int: 
        """
        __len__( arg1: VectSizeT) -> int

            C++ signature :
                unsigned long __len__(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: VectSizeT, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: VectSizeT, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: VectSizeT, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: VectSizeT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: VectSizeT, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: VectSizeT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<unsigned long, std::allocator<unsigned long> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class VectorOfStringVectors(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: VectorOfStringVectors, arg2: object) -> bool: 
        """
        __contains__( arg1: VectorOfStringVectors, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: VectorOfStringVectors, arg2: object) -> None: 
        """
        __delitem__( arg1: VectorOfStringVectors, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >*, std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > > > __iter__(boost::python::back_reference<std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > >&>)
        """
    @staticmethod
    def __len__( arg1: VectorOfStringVectors) -> int: 
        """
        __len__( arg1: VectorOfStringVectors) -> int

            C++ signature :
                unsigned long __len__(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: VectorOfStringVectors, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: VectorOfStringVectors, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: VectorOfStringVectors, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >, std::allocator<std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
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
def CreateDifferenceFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> UIntSparseIntVect:
    """
    CreateDifferenceFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> UIntSparseIntVect
        construct a difference fingerprint for a ChemicalReaction by subtracting the reactant fingerprint from the product fingerprint

        C++ signature :
            RDKit::SparseIntVect<unsigned int>* CreateDifferenceFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7f83635c8840>])
    """
def CreateStructuralFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> ExplicitBitVect:
    """
    CreateStructuralFingerprintForReaction( reaction: ChemicalReaction, ReactionFingerPrintParams: ReactionFingerprintParams = ReactionFingerprintParams()) -> ExplicitBitVect
        construct a structural fingerprint for a ChemicalReaction by concatenating the reactant fingerprint and the product fingerprint

        C++ signature :
            ExplicitBitVect* CreateStructuralFingerprintForReaction(RDKit::ChemicalReaction [,RDKit::ReactionFingerprintParams=<rdkit.Chem.rdChemReactions.ReactionFingerprintParams object at 0x7f83635c8240>])
    """
def EnumerateLibraryCanSerialize() -> bool:
    """
    EnumerateLibraryCanSerialize() -> bool
        Returns True if the EnumerateLibrary is serializable (requires boost serialization

        C++ signature :
            bool EnumerateLibraryCanSerialize()
    """
def GetChemDrawRxnAdjustParams() -> AdjustQueryParameters:
    """
    GetChemDrawRxnAdjustParams() -> AdjustQueryParameters
        (deprecated, see MatchOnlyAtRgroupsAdjustParams)
            Returns the chemdraw style adjustment parameters for reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters GetChemDrawRxnAdjustParams()
    """
def GetDefaultAdjustParams() -> AdjustQueryParameters:
    """
    GetDefaultAdjustParams() -> AdjustQueryParameters
        Returns the default adjustment parameters for reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters GetDefaultAdjustParams()
    """
def HasAgentTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    HasAgentTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool
        tests if the agents of a queryReaction are the same as those of a reaction

        C++ signature :
            bool HasAgentTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
    """
def HasProductTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool:
    """
    HasProductTemplateSubstructMatch( reaction: ChemicalReaction, queryReaction: ChemicalReaction) -> bool
        tests if the products of a queryReaction are substructures of the products of a reaction

        C++ signature :
            bool HasProductTemplateSubstructMatch(RDKit::ChemicalReaction,RDKit::ChemicalReaction)
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
def IsReactionTemplateMoleculeAgent( molecule: Mol, agentThreshold: float) -> bool:
    """
    IsReactionTemplateMoleculeAgent( molecule: Mol, agentThreshold: float) -> bool
        tests if a molecule can be classified as an agent depending on the ratio of mapped atoms and a give threshold

        C++ signature :
            bool IsReactionTemplateMoleculeAgent(RDKit::ROMol,double)
    """
def MatchOnlyAtRgroupsAdjustParams() -> AdjustQueryParameters:
    """
    MatchOnlyAtRgroupsAdjustParams() -> AdjustQueryParameters
        Only match at the specified rgroup locations in the reactant templates

        C++ signature :
            RDKit::MolOps::AdjustQueryParameters MatchOnlyAtRgroupsAdjustParams()
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
def ReactionToMrvFile( arg1: ChemicalReaction, reaction: str, prettyPrint: bool = False) -> None:
    """
    ReactionToMrvFile( arg1: ChemicalReaction, reaction: str, prettyPrint: bool = False) -> None
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
def ReduceProductToSideChains( product: Mol, addDummyAtoms: bool = True) -> Mol:
    """
    ReduceProductToSideChains( product: Mol, addDummyAtoms: bool = True) -> Mol
        reduce the product of a reaction to the side chains added by the reaction.              The output is a molecule with attached wildcards indicating where the product was attached.              The dummy atom has the same reaction-map number as the product atom (if available).

        C++ signature :
            RDKit::ROMol* ReduceProductToSideChains(boost::shared_ptr<RDKit::ROMol> [,bool=True])
    """
def RemoveMappingNumbersFromReactions( reaction: ChemicalReaction) -> None:
    """
    RemoveMappingNumbersFromReactions( reaction: ChemicalReaction) -> None
        Removes the mapping numbers from the molecules of a reaction

        C++ signature :
            void RemoveMappingNumbersFromReactions(RDKit::ChemicalReaction)
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
            RDKit::RxnOps::SanitizeRxnFlags SanitizeRxn(RDKit::ChemicalReaction {lvalue} [,unsigned long=4294967295 [,RDKit::MolOps::AdjustQueryParameters=<rdkit.Chem.rdmolops.AdjustQueryParameters object at 0x7f8363636360> [,bool=False]]])
    """
def UpdateProductsStereochemistry( reaction: ChemicalReaction) -> None:
    """
    UpdateProductsStereochemistry( reaction: ChemicalReaction) -> None
        Caution: This is an expert-user function which will change a property (molInversionFlag) of your products.          This function is called by default using the RXN or SMARTS parser for reactions and should really only be called if reactions have been constructed some other way.          The function updates the stereochemistry of the product by considering 4 different cases: inversion, retention, removal, and introduction

        C++ signature :
            void UpdateProductsStereochemistry(RDKit::ChemicalReaction*)
    """
SANITIZE_ADJUST_REACTANTS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
SANITIZE_ALL = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
SANITIZE_ATOM_MAPS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
SANITIZE_MERGEHS = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
SANITIZE_NONE = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
SANITIZE_RGROUP_NAMES = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
