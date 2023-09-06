"""Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm."""
from __future__ import annotations
import rdkit.Chem.rdRascalMCES
import typing
import Boost.Python

__all__ = [
    "FindMCES",
    "RascalButinaCluster",
    "RascalCluster",
    "RascalClusterOptions",
    "RascalOptions",
    "RascalResult"
]


class RascalClusterOptions(Boost.Python.instance):
    """
    RASCAL Cluster Options.  Most of these pertain to RascalCluster calculations.  Only similarityCutoff is used by RascalButinaCluster.
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def a(self) -> None:
        """
        The penalty score for each unconnected component in the MCES. Default=0.05.

        :type: None
        """
    @property
    def b(self) -> None:
        """
        The weight of matched bonds over matched atoms. Default=2.

        :type: None
        """
    @property
    def clusterMergeSim(self) -> None:
        """
        Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6.

        :type: None
        """
    @property
    def maxNumFrags(self) -> None:
        """
        The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity.

        :type: None
        """
    @property
    def minFragSize(self) -> None:
        """
        The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3.

        :type: None
        """
    @property
    def minIntraClusterSim(self) -> None:
        """
        Two pairs of molecules are included in the same cluster if the similarity between their MCESs is greater than this.  Default=0.9.

        :type: None
        """
    @property
    def numThreads(self) -> None:
        """
        Number of threads to use during clustering.  Default=-1 means all the hardware threads less one.

        :type: None
        """
    @property
    def similarityCutoff(self) -> None:
        """
        Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7.

        :type: None
        """
    __instance_size__ = 80
    pass
class RascalOptions(Boost.Python.instance):
    """
    RASCAL Options
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def allBestMCESs(self) -> None:
        """
        If True, reports all MCESs found of the same maximum size.  Default False means just report the first found.

        :type: None
        """
    @property
    def completeAromaticRings(self) -> None:
        """
        If True (default), partial aromatic rings won't be returned.

        :type: None
        """
    @property
    def maxFragSeparation(self) -> None:
        """
        Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed.

        :type: None
        """
    @property
    def minFragSize(self) -> None:
        """
        Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum.

        :type: None
        """
    @property
    def returnEmptyMCES(self) -> None:
        """
        If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned.

        :type: None
        """
    @property
    def ringMatchesRingOnly(self) -> None:
        """
        If True (default), ring bonds won't match ring bonds.

        :type: None
        """
    @property
    def similarityThreshold(self) -> None:
        """
        Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7.

        :type: None
        """
    @property
    def singleLargestFrag(self) -> None:
        """
        Return the just single largest fragment of the MCES.  This is equivalent to running with allBestMCEs=True, finding the result with the largest largestFragmentSize, and calling its largestFragmentOnly method.

        :type: None
        """
    @property
    def timeout(self) -> None:
        """
        Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit.

        :type: None
        """
    __instance_size__ = 64
    pass
class RascalResult(Boost.Python.instance):
    """
    Used to return RASCAL MCES results.
    """
    @staticmethod
    def atomMatches( arg1: RascalResult) -> list: 
        """
        atomMatches( arg1: RascalResult) -> list
            Likewise for atoms.

            C++ signature :
                boost::python::list atomMatches(RDKit::RascalMCES::RascalResult)
        """
    @staticmethod
    def bondMatches( arg1: RascalResult) -> list: 
        """
        bondMatches( arg1: RascalResult) -> list
            A function returning a list of list of tuples, each inner list containing the matching bonds in the MCES as tuples of bond indices from mol1 and mol2

            C++ signature :
                boost::python::list bondMatches(RDKit::RascalMCES::RascalResult)
        """
    @staticmethod
    def largestFragmentOnly( arg1: RascalResult) -> None: 
        """
        largestFragmentOnly( arg1: RascalResult) -> None
            Function that cuts the MCES down to the single largest frag.  This cannot be undone.

            C++ signature :
                void largestFragmentOnly(RDKit::RascalMCES::RascalResult {lvalue})
        """
    @property
    def largestFragmentSize(self) -> None:
        """
        Number of atoms in largest fragment.

        :type: None
        """
    @property
    def numFragments(self) -> None:
        """
        Number of fragments in MCES.

        :type: None
        """
    @property
    def similarity(self) -> None:
        """
        Johnson similarity between 2 molecules.

        :type: None
        """
    @property
    def smartsString(self) -> None:
        """
        SMARTS string defining the MCES.

        :type: None
        """
    @property
    def tier1Sim(self) -> None:
        """
        The tier 1 similarity estimate.

        :type: None
        """
    @property
    def tier2Sim(self) -> None:
        """
        The tier 2 similarity estimate.

        :type: None
        """
    @property
    def timedOut(self) -> None:
        """
        Whether it timed out.

        :type: None
        """
    pass
def FindMCES( mol1: Mol, mol2: Mol, opts: AtomPairsParameters = None) -> list:
    """
    FindMCES( mol1: Mol, mol2: Mol, opts: AtomPairsParameters = None) -> list
        Find one or more MCESs between the 2 molecules given.  Returns a list of RascalResult objects.- mol1- mol2 The two molecules for which to find the MCES- opts Optional RascalOptions object changing the default run mode.

        C++ signature :
            boost::python::list FindMCES(RDKit::ROMol,RDKit::ROMol [,boost::python::api::object=None])
    """
def RascalButinaCluster( mols: AtomPairsParameters, opts: AtomPairsParameters = None) -> list:
    """
    RascalButinaCluster( mols: AtomPairsParameters, opts: AtomPairsParameters = None) -> list
        Use the RASCAL MCES similarity metric to do Butina clustering (Butina JCICS 39 747-750 (1999)).  Returns a list of lists of molecules, each inner list being a cluster.  The last cluster is all the molecules that didn't fit into another cluster (the singletons).- mols List of molecules to be clustered- opts Optional RascalOptions object changing the default run mode.

        C++ signature :
            boost::python::list RascalButinaCluster(boost::python::api::object [,boost::python::api::object=None])
    """
def RascalCluster( mols: AtomPairsParameters, opts: AtomPairsParameters = None) -> list:
    """
    RascalCluster( mols: AtomPairsParameters, opts: AtomPairsParameters = None) -> list
        Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a list of lists of molecules, each inner list being a cluster.  The last cluster is all the molecules that didn't fit into another cluster (the singletons).- mols List of molecules to be clustered- opts Optional RascalOptions object changing the default run mode.

        C++ signature :
            boost::python::list RascalCluster(boost::python::api::object [,boost::python::api::object=None])
    """
