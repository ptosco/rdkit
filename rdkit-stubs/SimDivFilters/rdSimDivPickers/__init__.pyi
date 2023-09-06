"""Module containing the diversity and similarity pickers"""
from __future__ import annotations
import rdkit.SimDivFilters.rdSimDivPickers
import typing
import Boost.Python

__all__ = [
    "CENTROID",
    "CLINK",
    "ClusterMethod",
    "GOWER",
    "HierarchicalClusterPicker",
    "LeaderPicker",
    "MCQUITTY",
    "MaxMinPicker",
    "SLINK",
    "UPGMA",
    "WARD"
]


class ClusterMethod(Boost.Python.enum, int):
    CENTROID = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CENTROID
    CLINK = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CLINK
    GOWER = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.GOWER
    MCQUITTY = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.MCQUITTY
    SLINK = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.SLINK
    UPGMA = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.UPGMA
    WARD = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.WARD
    __slots__ = ()
    names = {'WARD': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.WARD, 'SLINK': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.SLINK, 'CLINK': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CLINK, 'UPGMA': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.UPGMA, 'MCQUITTY': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.MCQUITTY, 'GOWER': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.GOWER, 'CENTROID': rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CENTROID}
    values = {1: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.WARD, 2: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.SLINK, 3: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CLINK, 4: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.UPGMA, 5: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.MCQUITTY, 6: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.GOWER, 7: rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CENTROID}
    pass
class HierarchicalClusterPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using Hierarchical Clustering
    """
    @staticmethod
    def Cluster( arg1: HierarchicalClusterPicker, arg2: AtomPairsParameters, arg3: int, arg4: int) -> _vectSt6vectorIiSaIiEE: 
        """
        Cluster( arg1: HierarchicalClusterPicker, arg2: AtomPairsParameters, arg3: int, arg4: int) -> _vectSt6vectorIiSaIiEE
            Return a list of clusters of item from the pool using hierarchical clustering
            
            ARGUMENTS: 
              - distMat: 1D distance matrix (only the lower triangle elements)
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
            

            C++ signature :
                std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > Cluster(RDPickers::HierarchicalClusterPicker*,boost::python::api::object {lvalue},int,int)
        """
    @staticmethod
    def Pick( arg1: HierarchicalClusterPicker, arg2: AtomPairsParameters, arg3: int, arg4: int) -> _vecti: 
        """
        Pick( arg1: HierarchicalClusterPicker, arg2: AtomPairsParameters, arg3: int, arg4: int) -> _vecti
            Pick a diverse subset of items from a pool of items using hierarchical clustering
            
            ARGUMENTS: 
              - distMat: 1D distance matrix (only the lower triangle elements)
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
            

            C++ signature :
                std::vector<int, std::allocator<int> > Pick(RDPickers::HierarchicalClusterPicker*,boost::python::api::object {lvalue},int,int)
        """
    @staticmethod
    def __init__( arg1: object, clusterMethod: ClusterMethod) -> None: 
        """
        __init__( arg1: object, clusterMethod: ClusterMethod) -> None

            C++ signature :
                void __init__(_object*,RDPickers::HierarchicalClusterPicker::ClusterMethod)
        """
    __instance_size__ = 40
    pass
class LeaderPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using Roger Sayle's Leader algorithm (analogous to sphere exclusion). The algorithm is currently unpublished, but a description is available in this presentation from the 2019 RDKit UGM: https://github.com/rdkit/UGM_2019/raw/master/Presentations/Sayle_Clustering.pdf
    """
    def LazyBitVectorPick(self, objects: AtomPairsParameters, poolSize: int, threshold: float, pickSize: int = 0, firstPicks: AtomPairsParameters = (), numThreads: int = 1) -> _vecti: 
        """
        LazyBitVectorPick( self: LeaderPicker, objects: AtomPairsParameters, poolSize: int, threshold: float, pickSize: int = 0, firstPicks: AtomPairsParameters = (), numThreads: int = 1) -> _vecti
            Pick a subset of items from a collection of bit vectors using Tanimoto distance. The threshold value is a *distance* (i.e. 1-similarity). Note that the numThreads argument is currently ignored.

            C++ signature :
                std::vector<int, std::allocator<int> > LazyBitVectorPick(RDPickers::LeaderPicker*,boost::python::api::object,int,double [,int=0 [,boost::python::api::object=() [,int=1]]])
        """
    def LazyPick(self, distFunc: AtomPairsParameters, poolSize: int, threshold: float, pickSize: int = 0, firstPicks: AtomPairsParameters = (), numThreads: int = 1) -> _vecti: 
        """
        LazyPick( self: LeaderPicker, distFunc: AtomPairsParameters, poolSize: int, threshold: float, pickSize: int = 0, firstPicks: AtomPairsParameters = (), numThreads: int = 1) -> _vecti
            Pick a subset of items from a pool of items using the user-provided function to determine distances. Note that the numThreads argument is currently ignored.

            C++ signature :
                std::vector<int, std::allocator<int> > LazyPick(RDPickers::LeaderPicker*,boost::python::api::object,int,double [,int=0 [,boost::python::api::object=() [,int=1]]])
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 48
    pass
class MaxMinPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using the MaxMin Algorithm
    """
    def LazyBitVectorPick(self, objects: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1, useCache: AtomPairsParameters = None) -> _vecti: 
        """
        LazyBitVectorPick( self: MaxMinPicker, objects: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1, useCache: AtomPairsParameters = None) -> _vecti
            Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
            Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 
            ARGUMENTS:
            
              - vectors: a sequence of the bit vectors that should be picked from.
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
              - firstPicks: (optional) the first items to be picked (seeds the list)
              - seed: (optional) seed for the random number generator
              - useCache: IGNORED.
            

            C++ signature :
                std::vector<int, std::allocator<int> > LazyBitVectorPick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1 [,boost::python::api::object=None]]])
        """
    def LazyBitVectorPickWithThreshold(self, objects: AtomPairsParameters, poolSize: int, pickSize: int, threshold: float, firstPicks: AtomPairsParameters = (), seed: int = -1) -> tuple: 
        """
        LazyBitVectorPickWithThreshold( self: MaxMinPicker, objects: AtomPairsParameters, poolSize: int, pickSize: int, threshold: float, firstPicks: AtomPairsParameters = (), seed: int = -1) -> tuple
            Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
            Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 
            ARGUMENTS:
            
              - vectors: a sequence of the bit vectors that should be picked from.
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
              - threshold: stop picking when the distance goes below this value
              - firstPicks: (optional) the first items to be picked (seeds the list)
              - seed: (optional) seed for the random number generator
            

            C++ signature :
                boost::python::tuple LazyBitVectorPickWithThreshold(RDPickers::MaxMinPicker*,boost::python::api::object,int,int,double [,boost::python::api::object=() [,int=-1]])
        """
    def LazyPick(self, distFunc: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1, useCache: AtomPairsParameters = None) -> _vecti: 
        """
        LazyPick( self: MaxMinPicker, distFunc: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1, useCache: AtomPairsParameters = None) -> _vecti
            Pick a subset of items from a pool of items using the MaxMin Algorithm
            Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 
            ARGUMENTS:
            
              - distFunc: a function that should take two indices and return the
                          distance between those two points.
                          NOTE: the implementation caches distance values, so the
                          client code does not need to do so; indeed, it should not.
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
              - firstPicks: (optional) the first items to be picked (seeds the list)
              - seed: (optional) seed for the random number generator
              - useCache: IGNORED
            

            C++ signature :
                std::vector<int, std::allocator<int> > LazyPick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1 [,boost::python::api::object=None]]])
        """
    def LazyPickWithThreshold(self, distFunc: AtomPairsParameters, poolSize: int, pickSize: int, threshold: float, firstPicks: AtomPairsParameters = (), seed: int = -1) -> tuple: 
        """
        LazyPickWithThreshold( self: MaxMinPicker, distFunc: AtomPairsParameters, poolSize: int, pickSize: int, threshold: float, firstPicks: AtomPairsParameters = (), seed: int = -1) -> tuple
            Pick a subset of items from a pool of items using the MaxMin Algorithm
            Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 
            ARGUMENTS:
            
              - distFunc: a function that should take two indices and return the
                          distance between those two points.
                          NOTE: the implementation caches distance values, so the
                          client code does not need to do so; indeed, it should not.
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
              - threshold: stop picking when the distance goes below this value
              - firstPicks: (optional) the first items to be picked (seeds the list)
              - seed: (optional) seed for the random number generator
            

            C++ signature :
                boost::python::tuple LazyPickWithThreshold(RDPickers::MaxMinPicker*,boost::python::api::object,int,int,double [,boost::python::api::object=() [,int=-1]])
        """
    def Pick(self, distMat: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1) -> _vecti: 
        """
        Pick( self: MaxMinPicker, distMat: AtomPairsParameters, poolSize: int, pickSize: int, firstPicks: AtomPairsParameters = (), seed: int = -1) -> _vecti
            Pick a subset of items from a pool of items using the MaxMin Algorithm
            Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604 
            
            ARGUMENTS:
              - distMat: 1D distance matrix (only the lower triangle elements)
              - poolSize: number of items in the pool
              - pickSize: number of items to pick from the pool
              - firstPicks: (optional) the first items to be picked (seeds the list)
              - seed: (optional) seed for the random number generator
            

            C++ signature :
                std::vector<int, std::allocator<int> > Pick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1]])
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
CENTROID = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CENTROID
CLINK = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.CLINK
GOWER = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.GOWER
MCQUITTY = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.MCQUITTY
SLINK = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.SLINK
UPGMA = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.UPGMA
WARD = rdkit.SimDivFilters.rdSimDivPickers.ClusterMethod.WARD
