"""Module containing bunch of functions for information metrics and a ranker to rank bits"""
from __future__ import annotations
import rdkit.ML.InfoTheory.rdInfoTheory
import typing
import Boost.Python

__all__ = [
    "BIASCHISQUARE",
    "BIASENTROPY",
    "BitCorrMatGenerator",
    "CHISQUARE",
    "ChiSquare",
    "ENTROPY",
    "InfoBitRanker",
    "InfoEntropy",
    "InfoGain",
    "InfoType"
]


class BitCorrMatGenerator(Boost.Python.instance):
    """
    A class to generate a pairwise correlation matrix between a list of bits
    The mode of operation for this class is something like this

       >>> cmg = BitCorrMatGenerator() 
       >>> cmg.SetBitList(blist) 
       >>> for fp in fpList:  
       >>>    cmg.CollectVotes(fp)  
       >>> corrMat = cmg.GetCorrMatrix() 
        
       The resulting correlation matrix is a one dimensional nummeric array containing the 
       lower triangle elements
    """
    @staticmethod
    def CollectVotes( arg1: BitCorrMatGenerator, arg2: AtomPairsParameters) -> None: 
        """
        CollectVotes( arg1: BitCorrMatGenerator, arg2: AtomPairsParameters) -> None
            For each pair of on bits (bi, bj) in fp increase the correlation count for the pair by 1
            
            ARGUMENTS:
            
              - fp : a bit vector to collect the fingerprints from
            

            C++ signature :
                void CollectVotes(RDInfoTheory::BitCorrMatGenerator*,boost::python::api::object)
        """
    @staticmethod
    def GetCorrMatrix( arg1: BitCorrMatGenerator) -> object: 
        """
        GetCorrMatrix( arg1: BitCorrMatGenerator) -> object
            Get the correlation matrix following the collection of votes from a bunch of fingerprints
            

            C++ signature :
                _object* GetCorrMatrix(RDInfoTheory::BitCorrMatGenerator*)
        """
    @staticmethod
    def SetBitList( arg1: BitCorrMatGenerator, arg2: AtomPairsParameters) -> None: 
        """
        SetBitList( arg1: BitCorrMatGenerator, arg2: AtomPairsParameters) -> None
            Set the list of bits that need to be correllated
            
             This may for example be their top ranking ensemble bits
            
            ARGUMENTS:
            
              - bitList : an integer list of bit IDs
            

            C++ signature :
                void SetBitList(RDInfoTheory::BitCorrMatGenerator*,boost::python::api::object)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 64
    pass
class InfoBitRanker(Boost.Python.instance):
    """
    A class to rank the bits from a series of labelled fingerprints
    A simple demonstration may help clarify what this class does. 
    Here's a small set of vectors:

    >>> for i,bv in enumerate(bvs): print(bv.ToBitString(),acts[i])
    ... 
    0001 0
    0101 0
    0010 1
    1110 1

    Default ranker, using infogain:

    >>> ranker = InfoBitRanker(4,2)  
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ... 
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ... 
    3 1.000 2 0
    2 1.000 0 2
    0 0.311 0 1

    Using the biased infogain:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASENTROPY)
    >>> ranker.SetBiasList((1,))
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ... 
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ... 
    2 1.000 0 2
    0 0.311 0 1
    1 0.000 1 1

    A chi squared ranker is also available:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.CHISQUARE)
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ... 
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ... 
    3 4.000 2 0
    2 4.000 0 2
    0 1.333 0 1

    As is a biased chi squared:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASCHISQUARE)
    >>> ranker.SetBiasList((1,))
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ... 
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ... 
    2 4.000 0 2
    0 1.333 0 1
    1 0.000 1 1
    """
    @staticmethod
    def AccumulateVotes( arg1: InfoBitRanker, arg2: AtomPairsParameters, arg3: int) -> None: 
        """
        AccumulateVotes( arg1: InfoBitRanker, arg2: AtomPairsParameters, arg3: int) -> None
            Accumulate the votes for all the bits turned on in a bit vector
            
            ARGUMENTS:
            
              - bv : bit vector either ExplicitBitVect or SparseBitVect operator
              - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses 
            

            C++ signature :
                void AccumulateVotes(RDInfoTheory::InfoBitRanker*,boost::python::api::object,int)
        """
    @staticmethod
    def GetTopN( arg1: InfoBitRanker, arg2: int) -> object: 
        """
        GetTopN( arg1: InfoBitRanker, arg2: int) -> object
            Returns the top n bits ranked by the information metric
            This is actually the function where most of the work of ranking is happening
            
            ARGUMENTS:
            
              - num : the number of top ranked bits that are required
            

            C++ signature :
                _object* GetTopN(RDInfoTheory::InfoBitRanker*,int)
        """
    @staticmethod
    def SetBiasList( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None: 
        """
        SetBiasList( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None
            Set the classes to which the entropy calculation should be biased
            
            This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits. 
            In this mode, a bit must be correlated higher with one of the biased classes than all the 
            other classes. For example, in a two class problem with actives and inactives, the fraction of 
            actives that hit the bit has to be greater than the fraction of inactives that hit the bit
            
            ARGUMENTS: 
            
              - classList : list of class ids that we want a bias towards
            

            C++ signature :
                void SetBiasList(RDInfoTheory::InfoBitRanker*,boost::python::api::object)
        """
    @staticmethod
    def SetMaskBits( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None: 
        """
        SetMaskBits( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None
            Set the mask bits for the calculation
            
            ARGUMENTS: 
            
              - maskBits : list of mask bits to use
            

            C++ signature :
                void SetMaskBits(RDInfoTheory::InfoBitRanker*,boost::python::api::object)
        """
    @staticmethod
    def Tester( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None: 
        """
        Tester( arg1: InfoBitRanker, arg2: AtomPairsParameters) -> None

            C++ signature :
                void Tester(RDInfoTheory::InfoBitRanker*,boost::python::api::object)
        """
    @staticmethod
    def WriteTopBitsToFile( arg1: InfoBitRanker, arg2: str) -> None: 
        """
        WriteTopBitsToFile( arg1: InfoBitRanker, arg2: str) -> None
            Write the bits that have been ranked to a file

            C++ signature :
                void WriteTopBitsToFile(RDInfoTheory::InfoBitRanker {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, nBits: int, nClasses: int) -> None: 
        """
        __init__( arg1: object, nBits: int, nClasses: int) -> None

            C++ signature :
                void __init__(_object*,int,int)

            C++ signature :
                void __init__(_object*,int,int,RDInfoTheory::InfoBitRanker::InfoType)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, nBits: int, nClasses: int, infoType: InfoType) -> None: ...
    __instance_size__ = 136
    pass
class InfoType(Boost.Python.enum, int):
    BIASCHISQUARE = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASCHISQUARE
    BIASENTROPY = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASENTROPY
    CHISQUARE = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.CHISQUARE
    ENTROPY = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.ENTROPY
    __slots__ = ()
    names = {'ENTROPY': rdkit.ML.InfoTheory.rdInfoTheory.InfoType.ENTROPY, 'BIASENTROPY': rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASENTROPY, 'CHISQUARE': rdkit.ML.InfoTheory.rdInfoTheory.InfoType.CHISQUARE, 'BIASCHISQUARE': rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASCHISQUARE}
    values = {1: rdkit.ML.InfoTheory.rdInfoTheory.InfoType.ENTROPY, 2: rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASENTROPY, 3: rdkit.ML.InfoTheory.rdInfoTheory.InfoType.CHISQUARE, 4: rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASCHISQUARE}
    pass
def ChiSquare( arg1: AtomPairsParameters) -> float:
    """
    ChiSquare( arg1: AtomPairsParameters) -> float
        Calculates the chi squared value for a variable
        
           ARGUMENTS:
        
             - varMat: a Numeric Array object
               varMat is a Numeric array with the number of possible occurrences
                 of each result for reach possible value of the given variable.
        
               So, for a variable which adopts 4 possible values and a result which
                 has 3 possible values, varMat would be 4x3
        
           RETURNS:
        
             - a Python float object
        

        C++ signature :
            double ChiSquare(boost::python::api::object)
    """
def InfoEntropy( arg1: AtomPairsParameters) -> float:
    """
    InfoEntropy( arg1: AtomPairsParameters) -> float
        calculates the informational entropy of the values in an array
        
          ARGUMENTS:
            
            - resMat: pointer to a long int array containing the data
            - dim: long int containing the length of the _tPtr_ array.
        
          RETURNS:
        
            a double
        

        C++ signature :
            double InfoEntropy(boost::python::api::object)
    """
def InfoGain( arg1: AtomPairsParameters) -> float:
    """
    InfoGain( arg1: AtomPairsParameters) -> float
        Calculates the information gain for a variable
        
           ARGUMENTS:
        
             - varMat: a Numeric Array object
               varMat is a Numeric array with the number of possible occurrences
                 of each result for reach possible value of the given variable.
        
               So, for a variable which adopts 4 possible values and a result which
                 has 3 possible values, varMat would be 4x3
        
           RETURNS:
        
             - a Python float object
        
           NOTES
        
             - this is a dropin replacement for _PyInfoGain()_ in entropy.py
        

        C++ signature :
            double InfoGain(boost::python::api::object)
    """
BIASCHISQUARE = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASCHISQUARE
BIASENTROPY = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASENTROPY
CHISQUARE = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.CHISQUARE
ENTROPY = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.ENTROPY
