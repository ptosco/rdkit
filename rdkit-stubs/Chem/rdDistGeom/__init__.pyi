"""Module containing functions to compute atomic coordinates in 3D using distance geometry"""
from __future__ import annotations
import rdkit.Chem.rdDistGeom
import typing
import Boost.Python

__all__ = [
    "BAD_DOUBLE_BOND_STEREO",
    "CHECK_CHIRAL_CENTERS",
    "CHECK_TETRAHEDRAL_CENTERS",
    "ETDG",
    "ETKDG",
    "ETKDGv2",
    "ETKDGv3",
    "ETK_MINIMIZATION",
    "EmbedFailureCauses",
    "EmbedMolecule",
    "EmbedMultipleConfs",
    "EmbedParameters",
    "FINAL_CENTER_IN_VOLUME",
    "FINAL_CHIRAL_BOUNDS",
    "FIRST_MINIMIZATION",
    "GetExperimentalTorsions",
    "GetMoleculeBoundsMatrix",
    "INITIAL_COORDS",
    "KDG",
    "LINEAR_DOUBLE_BOND",
    "MINIMIZE_FOURTH_DIMENSION",
    "srETKDGv3"
]


class EmbedFailureCauses(Boost.Python.enum, int):
    BAD_DOUBLE_BOND_STEREO = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
    CHECK_CHIRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
    CHECK_TETRAHEDRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
    ETK_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
    FINAL_CENTER_IN_VOLUME = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
    FINAL_CHIRAL_BOUNDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
    FIRST_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
    INITIAL_COORDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
    LINEAR_DOUBLE_BOND = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
    MINIMIZE_FOURTH_DIMENSION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
    __slots__ = ()
    names = {'INITIAL_COORDS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS, 'FIRST_MINIMIZATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION, 'CHECK_TETRAHEDRAL_CENTERS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS, 'CHECK_CHIRAL_CENTERS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS, 'MINIMIZE_FOURTH_DIMENSION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION, 'ETK_MINIMIZATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION, 'FINAL_CHIRAL_BOUNDS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS, 'FINAL_CENTER_IN_VOLUME': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME, 'LINEAR_DOUBLE_BOND': rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND, 'BAD_DOUBLE_BOND_STEREO': rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO}
    values = {0: rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS, 1: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION, 2: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS, 3: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS, 4: rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION, 5: rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION, 6: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS, 7: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME, 8: rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND, 9: rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO}
    pass
class EmbedParameters(Boost.Python.instance):
    """
    Parameters controlling embedding
    """
    @staticmethod
    def GetFailureCounts( arg1: EmbedParameters) -> tuple: 
        """
        GetFailureCounts( arg1: EmbedParameters) -> tuple
            returns the counts of each failure type

            C++ signature :
                boost::python::tuple GetFailureCounts(RDKit::DGeomHelpers::EmbedParameters*)
        """
    @staticmethod
    def SetBoundsMat( arg1: EmbedParameters, arg2: AtomPairsParameters) -> None: 
        """
        SetBoundsMat( arg1: EmbedParameters, arg2: AtomPairsParameters) -> None
            set the distance-bounds matrix to be used (no triangle smoothing will be done on this) from a Numpy array

            C++ signature :
                void SetBoundsMat(RDKit::DGeomHelpers::EmbedParameters*,boost::python::api::object)
        """
    @staticmethod
    def SetCPCI( arg1: EmbedParameters, arg2: dict) -> None: 
        """
        SetCPCI( arg1: EmbedParameters, arg2: dict) -> None
            set the customised pairwise Columb-like interaction to atom pairs.used during structural minimisation stage

            C++ signature :
                void SetCPCI(RDKit::DGeomHelpers::EmbedParameters*,boost::python::dict {lvalue})
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def ETversion(self) -> None:
        """
        version of the experimental torsion-angle preferences

        :type: None
        """
    @property
    def boundsMatForceScaling(self) -> None:
        """
        scale the weights of the atom pair distance restraints relative to the other types of restraints

        :type: None
        """
    @property
    def boxSizeMult(self) -> None:
        """
        determines the size of the box used for random coordinates

        :type: None
        """
    @property
    def clearConfs(self) -> None:
        """
        clear all existing conformations on the molecule

        :type: None
        """
    @property
    def embedFragmentsSeparately(self) -> None:
        """
        split the molecule into fragments and embed them separately

        :type: None
        """
    @property
    def enableSequentialRandomSeeds(self) -> None:
        """
        handle random number seeds so that conformer generation can be restarted

        :type: None
        """
    @property
    def enforceChirality(self) -> None:
        """
        enforce correct chirilaty if chiral centers are present

        :type: None
        """
    @property
    def forceTransAmides(self) -> None:
        """
        constrain amide bonds to be trans

        :type: None
        """
    @property
    def ignoreSmoothingFailures(self) -> None:
        """
        try and embed the molecule if if triangle smoothing of the bounds matrix fails

        :type: None
        """
    @property
    def maxIterations(self) -> None:
        """
        maximum number of embedding attempts to use for a single conformation

        :type: None
        """
    @property
    def numThreads(self) -> None:
        """
        number of threads to use when embedding multiple conformations

        :type: None
        """
    @property
    def numZeroFail(self) -> None:
        """
        fail embedding if we have at least this many zero eigenvalues

        :type: None
        """
    @property
    def onlyHeavyAtomsForRMS(self) -> None:
        """
        Only consider heavy atoms when doing RMS filtering

        :type: None
        """
    @property
    def optimizerForceTol(self) -> None:
        """
        the tolerance to be used during the distance-geometry force field minimization

        :type: None
        """
    @property
    def pruneRmsThresh(self) -> None:
        """
        used to filter multiple conformations: keep only conformations that are at least this far apart from each other

        :type: None
        """
    @property
    def randNegEig(self) -> None:
        """
        if the embedding yields a negative eigenvalue, pick coordinates that correspond to this component at random

        :type: None
        """
    @property
    def randomSeed(self) -> None:
        """
        seed for the random number generator

        :type: None
        """
    @property
    def trackFailures(self) -> None:
        """
        keep track of which checks during the embedding process fail

        :type: None
        """
    @property
    def useBasicKnowledge(self) -> None:
        """
        impose basic-knowledge constraints such as flat rings

        :type: None
        """
    @property
    def useExpTorsionAnglePrefs(self) -> None:
        """
        impose experimental torsion angle preferences

        :type: None
        """
    @property
    def useMacrocycleTorsions(self) -> None:
        """
        impose macrocycle torsion angle preferences

        :type: None
        """
    @property
    def useRandomCoords(self) -> None:
        """
        start the embedding from random coordinates instead of using eigenvalues of the distance matrix

        :type: None
        """
    @property
    def useSmallRingTorsions(self) -> None:
        """
        impose small ring torsion angle preferences

        :type: None
        """
    @property
    def useSymmetryForPruning(self) -> None:
        """
        use molecule symmetry when doing the RMSD pruning. Note that this option automatically also sets onlyHeavyAtomsForRMS to true.

        :type: None
        """
    @property
    def verbose(self) -> None:
        """
        be verbose about configuration

        :type: None
        """
    __instance_size__ = 208
    pass
def ETDG() -> EmbedParameters:
    """
    ETDG() -> EmbedParameters
        Returns an EmbedParameters object for the ETDG method.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETDG()
    """
def ETKDG() -> EmbedParameters:
    """
    ETKDG() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 1.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDG()
    """
def ETKDGv2() -> EmbedParameters:
    """
    ETKDGv2() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 2.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDGv2()
    """
def ETKDGv3() -> EmbedParameters:
    """
    ETKDGv3() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 3 (macrocycles).

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* ETKDGv3()
    """
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
@typing.overload
def EmbedMultipleConfs( mol: Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> _vecti:
    """
    EmbedMultipleConfs( mol: Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = False, ETversion: int = 1) -> _vecti
        Use distance geometry to obtain multiple sets of 
         coordinates for a molecule
         
         ARGUMENTS:
        
          - mol : the molecule of interest
          - numConfs : the number of conformers to generate 
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
          - pruneRmsThresh : Retain only the conformations out of 'numConfs' 
                            after embedding that are at least 
                            this far apart from each other. 
                            RMSD is computed on the heavy atoms. 
                            Pruning is greedy; i.e. the first embedded conformation
                            is retained and from then on only those that are at
                            least pruneRmsThresh away from all retained conformations
                            are kept. The pruning is done after embedding and 
                            bounds violation minimization. No pruning by default.
          - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                       require some atoms to have fixed coordinates in the resulting 
                       conformation.
          - forceTol : tolerance to be used during the force-field minimization with 
                       the distance geometry force field.
          - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                       of the bounds matrix fails.
          - enforceChirality : enforce the correct chirality if chiral centers are present.
          - numThreads : number of threads to use while embedding. This only has an effect if the RDKit
                       was built with multi-thread support.
                      If set to zero, the max supported by the system will be used.
          - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
          - useBasicKnowledge : impose basic knowledge such as flat rings
          - printExpTorsionAngles : print the output from the experimental torsion angles
         RETURNS:
        
            List of new conformation IDs 
        
        

        C++ signature :
            std::vector<int, std::allocator<int> > EmbedMultipleConfs(RDKit::ROMol {lvalue} [,unsigned int=10 [,unsigned int=0 [,int=-1 [,bool=True [,bool=False [,double=2.0 [,bool=True [,unsigned int=1 [,double=-1.0 [,boost::python::dict {lvalue}={} [,double=0.001 [,bool=False [,bool=True [,int=1 [,bool=True [,bool=True [,bool=False [,bool=False [,bool=False [,unsigned int=1]]]]]]]]]]]]]]]]]]]])

        C++ signature :
            std::vector<int, std::allocator<int> > EmbedMultipleConfs(RDKit::ROMol {lvalue},unsigned int,RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
@typing.overload
def EmbedMultipleConfs( mol: Mol, numConfs: int, params: EmbedParameters) -> _vecti:
    pass
@typing.overload
def GetExperimentalTorsions( mol: Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> tuple:
    """
    GetExperimentalTorsions( mol: Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> tuple
        returns information about the bonds corresponding to experimental torsions

        C++ signature :
            boost::python::tuple GetExperimentalTorsions(RDKit::ROMol [,bool=True [,bool=False [,bool=True [,bool=True [,unsigned int=2 [,bool=False]]]]]])

        C++ signature :
            boost::python::tuple GetExperimentalTorsions(RDKit::ROMol,RDKit::DGeomHelpers::EmbedParameters)
    """
@typing.overload
def GetExperimentalTorsions( mol: Mol, embedParams: EmbedParameters) -> tuple:
    pass
def GetMoleculeBoundsMatrix( mol: Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False) -> object:
    """
    GetMoleculeBoundsMatrix( mol: Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False) -> object
        Returns the distance bounds matrix for a molecule
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - set15bounds : set bounds for 1-5 atom distances based on 
                            topology (otherwise stop at 1-4s)
            - scaleVDW : scale down the sum of VDW radii when setting the 
                         lower bounds for atoms less than 5 bonds apart 
            - doTriangleSmoothing : run triangle smoothing on the bounds 
                         matrix before returning it 
         RETURNS:
        
            the bounds matrix as a Numeric array with lower bounds in 
            the lower triangle and upper bounds in the upper triangle
        
        

        C++ signature :
            _object* GetMoleculeBoundsMatrix(RDKit::ROMol {lvalue} [,bool=True [,bool=False [,bool=True [,bool=False]]]])
    """
def KDG() -> EmbedParameters:
    """
    KDG() -> EmbedParameters
        Returns an EmbedParameters object for the KDG method.

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* KDG()
    """
def srETKDGv3() -> EmbedParameters:
    """
    srETKDGv3() -> EmbedParameters
        Returns an EmbedParameters object for the ETKDG method - version 3 (small rings).

        C++ signature :
            RDKit::DGeomHelpers::EmbedParameters* srETKDGv3()
    """
BAD_DOUBLE_BOND_STEREO = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
CHECK_CHIRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
CHECK_TETRAHEDRAL_CENTERS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
ETK_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
FINAL_CENTER_IN_VOLUME = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
FINAL_CHIRAL_BOUNDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
FIRST_MINIMIZATION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
INITIAL_COORDS = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
LINEAR_DOUBLE_BOND = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
MINIMIZE_FOURTH_DIMENSION = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
