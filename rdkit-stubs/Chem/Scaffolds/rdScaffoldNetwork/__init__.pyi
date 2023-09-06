"""Module containing functions for creating a Scaffold Network"""
from __future__ import annotations
import rdkit.Chem.Scaffolds.rdScaffoldNetwork
import typing
import Boost.Python

__all__ = [
    "BRICSScaffoldParams",
    "CreateScaffoldNetwork",
    "EdgeType",
    "NetworkEdge",
    "NetworkEdge_VECT",
    "ScaffoldNetwork",
    "ScaffoldNetworkParams",
    "UpdateScaffoldNetwork"
]


class EdgeType(Boost.Python.enum, int):
    Fragment = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment
    Generic = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic
    GenericBond = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond
    Initialize = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize
    RemoveAttachment = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment
    __slots__ = ()
    names = {'Fragment': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment, 'Generic': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic, 'GenericBond': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond, 'RemoveAttachment': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment, 'Initialize': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize}
    values = {1: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment, 2: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic, 3: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond, 4: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment, 5: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize}
    pass
class NetworkEdge(Boost.Python.instance):
    """
    A scaffold network edge
    """
    @staticmethod
    def __str__( arg1: NetworkEdge) -> object: 
        """
        __str__( arg1: NetworkEdge) -> object

            C++ signature :
                _object* __str__(RDKit::ScaffoldNetwork::NetworkEdge {lvalue})
        """
    @property
    def beginIdx(self) -> None:
        """
        index of the begin node in node list

        :type: None
        """
    @property
    def endIdx(self) -> None:
        """
        index of the end node in node list

        :type: None
        """
    @property
    def type(self) -> None:
        """
        type of the edge

        :type: None
        """
    pass
class NetworkEdge_VECT(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: NetworkEdge_VECT, arg2: object) -> bool: 
        """
        __contains__( arg1: NetworkEdge_VECT, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: NetworkEdge_VECT, arg2: object) -> None: 
        """
        __delitem__( arg1: NetworkEdge_VECT, arg2: object) -> None

            C++ signature :
                void __delitem__(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> >&>,_object*)
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, __gnu_cxx::__normal_iterator<RDKit::ScaffoldNetwork::NetworkEdge*, std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > > > __iter__(boost::python::back_reference<std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> >&>)
        """
    @staticmethod
    def __len__( arg1: NetworkEdge_VECT) -> int: 
        """
        __len__( arg1: NetworkEdge_VECT) -> int

            C++ signature :
                unsigned long __len__(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: NetworkEdge_VECT, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: NetworkEdge_VECT, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},_object*,_object*)
        """
    @staticmethod
    def append( arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None: 
        """
        append( arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void append(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},boost::python::api::object)
        """
    @staticmethod
    def extend( arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None: 
        """
        extend( arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None

            C++ signature :
                void extend(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},boost::python::api::object)
        """
    __instance_size__ = 48
    pass
class ScaffoldNetwork(Boost.Python.instance):
    """
    A scaffold network
    """
    @staticmethod
    def __getinitargs__( arg1: ScaffoldNetwork) -> tuple: 
        """
        __getinitargs__( arg1: ScaffoldNetwork) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDKit::ScaffoldNetwork::ScaffoldNetwork)
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
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @property
    def counts(self) -> None:
        """
        the number of times each node was encountered while building the network.

        :type: None
        """
    @property
    def edges(self) -> None:
        """
        the sequence of network edges

        :type: None
        """
    @property
    def molCounts(self) -> None:
        """
        the number of moleclues each node was found in.

        :type: None
        """
    @property
    def nodes(self) -> None:
        """
        the sequence of SMILES defining the nodes

        :type: None
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 120
    __safe_for_unpickling__ = True
    pass
class ScaffoldNetworkParams(Boost.Python.instance):
    """
    Scaffold network parameters
    """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, bondBreakerSmartsList: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE) -> None: ...
    @property
    def collectMolCounts(self) -> None:
        """
        keep track of the number of molecules each scaffold was found in

        :type: None
        """
    @property
    def flattenChirality(self) -> None:
        """
        remove chirality and bond stereo when flattening

        :type: None
        """
    @property
    def flattenIsotopes(self) -> None:
        """
        remove isotopes when flattening

        :type: None
        """
    @property
    def flattenKeepLargest(self) -> None:
        """
        keep only the largest fragment when doing flattening

        :type: None
        """
    @property
    def includeGenericBondScaffolds(self) -> None:
        """
        include scaffolds with all bonds replaced by single bonds

        :type: None
        """
    @property
    def includeGenericScaffolds(self) -> None:
        """
        include scaffolds with all atoms replaced by dummies

        :type: None
        """
    @property
    def includeScaffoldsWithAttachments(self) -> None:
        """
        Include the version of the scaffold with attachment points

        :type: None
        """
    @property
    def includeScaffoldsWithoutAttachments(self) -> None:
        """
        remove attachment points from scaffolds and include the result

        :type: None
        """
    @property
    def keepOnlyFirstFragment(self) -> None:
        """
        keep only the first fragment from the bond breaking rule

        :type: None
        """
    @property
    def pruneBeforeFragmenting(self) -> None:
        """
        Do a pruning/flattening step before starting fragmenting

        :type: None
        """
    __instance_size__ = 64
    pass
def BRICSScaffoldParams() -> ScaffoldNetworkParams:
    """
    BRICSScaffoldParams() -> ScaffoldNetworkParams
        Returns parameters for generating scaffolds using BRICS fragmentation rules

        C++ signature :
            RDKit::ScaffoldNetwork::ScaffoldNetworkParams* BRICSScaffoldParams()
    """
def CreateScaffoldNetwork( mols: AtomPairsParameters, params: ScaffoldNetworkParams) -> ScaffoldNetwork:
    """
    CreateScaffoldNetwork( mols: AtomPairsParameters, params: ScaffoldNetworkParams) -> ScaffoldNetwork
        create (and return) a new network from a sequence of molecules

        C++ signature :
            RDKit::ScaffoldNetwork::ScaffoldNetwork* CreateScaffoldNetwork(boost::python::api::object,RDKit::ScaffoldNetwork::ScaffoldNetworkParams)
    """
def UpdateScaffoldNetwork( mols: AtomPairsParameters, network: ScaffoldNetwork, params: ScaffoldNetworkParams) -> None:
    """
    UpdateScaffoldNetwork( mols: AtomPairsParameters, network: ScaffoldNetwork, params: ScaffoldNetworkParams) -> None
        update an existing network by adding molecules

        C++ signature :
            void UpdateScaffoldNetwork(boost::python::api::object,RDKit::ScaffoldNetwork::ScaffoldNetwork {lvalue},RDKit::ScaffoldNetwork::ScaffoldNetworkParams)
    """
