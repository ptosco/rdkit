from __future__ import annotations
import rdkit.Chem.EnumerateStereoisomers
import typing
import random
import rdkit.Chem

__all__ = [
    "Chem",
    "EmbedMolecule",
    "EnumerateStereoisomers",
    "GetStereoisomerCount",
    "StereoEnumerationOptions",
    "random"
]


class StereoEnumerationOptions():
    """
    - tryEmbedding: if set the process attempts to generate a standard RDKit distance geometry
      conformation for the stereisomer. If this fails, we assume that the stereoisomer is
      non-physical and don't return it. NOTE that this is computationally expensive and is
      just a heuristic that could result in stereoisomers being lost.

    - onlyUnassigned: if set (the default), stereocenters which have specified stereochemistry
      will not be perturbed unless they are part of a relative stereo
      group.

    - maxIsomers: the maximum number of isomers to yield, if the
      number of possible isomers is greater than maxIsomers, a
      random subset will be yielded. If 0, all isomers are
      yielded. Since every additional stereo center doubles the
      number of results (and execution time) it's important to
      keep an eye on this.

    - onlyStereoGroups: Only find stereoisomers that differ at the
      StereoGroups associated with the molecule.
    """
    __slots__ = ('tryEmbedding', 'onlyUnassigned', 'onlyStereoGroups', 'maxIsomers', 'rand', 'unique')
    maxIsomers: member_descriptor # value = <member 'maxIsomers' of 'StereoEnumerationOptions' objects>
    onlyStereoGroups: member_descriptor # value = <member 'onlyStereoGroups' of 'StereoEnumerationOptions' objects>
    onlyUnassigned: member_descriptor # value = <member 'onlyUnassigned' of 'StereoEnumerationOptions' objects>
    rand: member_descriptor # value = <member 'rand' of 'StereoEnumerationOptions' objects>
    tryEmbedding: member_descriptor # value = <member 'tryEmbedding' of 'StereoEnumerationOptions' objects>
    unique: member_descriptor # value = <member 'unique' of 'StereoEnumerationOptions' objects>
    pass
class _AtomFlipper():
    pass
class _BondFlipper():
    pass
class _RangeBitsGenerator():
    pass
class _StereoGroupFlipper():
    pass
class _UniqueRandomBitsGenerator():
    pass
@typing.overload
def EmbedMolecule( mol: Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> int:
    """
    EmbedMolecule( mol: Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> int
        Use distance geometry to obtain initial 
         coordinates for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - maxAttempts : the maximum number of attempts to try embedding 
            - randomSeed : provide a seed for the random number generator 
                           so that the same coordinates can be obtained 
                           for a molecule on multiple runs. If -1, the 
                           RNG will not be seeded. 
            - clearConfs : clear all existing conformations on the molecule
            - useRandomCoords : Start the embedding from random coordinates instead of
                                using eigenvalues of the distance matrix.
            - boxSizeMult    Determines the size of the box that is used for
                             random coordinates. If this is a positive number, the 
                             side length will equal the largest element of the distance
                             matrix times boxSizeMult. If this is a negative number,
                             the side length will equal -boxSizeMult (i.e. independent
                             of the elements of the distance matrix).
            - randNegEig : If the embedding yields a negative eigenvalue, 
                           pick coordinates that correspond 
                           to this component at random 
            - numZeroFail : fail embedding if we have at least this many zero eigenvalues 
            - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                         require some atoms to have fixed coordinates in the resulting 
                         conformation.
            - forceTol : tolerance to be used during the force-field minimization with 
                         the distance geometry force field.
            - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                         of the bounds matrix fails.
            - enforceChirality : enforce the correct chirality if chiral centers are present.
            - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
            - useBasicKnowledge : impose basic knowledge such as flat rings
            - printExpTorsionAngles : print the output from the experimental torsion angles
        
         RETURNS:
        
            ID of the new conformation added to the molecule 
        
        

        C++ signature :
            int EmbedMolecule(RDKit::ROMol {lvalue} [,unsigned int=0 [,int=-1 [,bool=True [,bool=False [,double=2.0 [,bool=True [,unsigned int=1 [,boost::python::dict {lvalue}={} [,double=0.001 [,bool=False [,bool=True [,bool=True [,bool=True [,bool=False [,bool=False [,bool=False [,unsigned int=1]]]]]]]]]]]]]]]]])

        C++ signature :
            int EmbedMolecule(RDKit::ROMol {lvalue},RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
@typing.overload
def EmbedMolecule( mol: Mol, params: EmbedParameters) -> int:
    pass
