"""Module containing Trajectory and Snapshot objects
"""
from __future__ import annotations
import rdkit.Chem.rdtrajectory
import typing
import Boost.Python

__all__ = [
    "ReadAmberTrajectory",
    "ReadGromosTrajectory",
    "Snapshot",
    "Trajectory"
]


class Snapshot(Boost.Python.instance):
    """
    A class which allows storing coordinates from a trajectory
    """
    def GetEnergy(self) -> float: 
        """
        GetEnergy( self: Snapshot) -> float
            returns the energy for this Snapshot

            C++ signature :
                double GetEnergy(RDKit::Snapshot {lvalue})
        """
    def GetPoint2D(self, pointNum: int) -> Point2D: 
        """
        GetPoint2D( self: Snapshot, pointNum: int) -> Point2D
            return the coordinates at pointNum as a Point2D object; requires the Trajectory dimension to be == 2

            C++ signature :
                RDGeom::Point2D GetPoint2D(RDKit::Snapshot {lvalue},unsigned int)
        """
    def GetPoint3D(self, pointNum: int) -> Point3D: 
        """
        GetPoint3D( self: Snapshot, pointNum: int) -> Point3D
            return the coordinates at pointNum as a Point3D object; requires the Trajectory dimension to be >= 2

            C++ signature :
                RDGeom::Point3D GetPoint3D(RDKit::Snapshot {lvalue},unsigned int)
        """
    def SetEnergy(self, energy: float) -> None: 
        """
        SetEnergy( self: Snapshot, energy: float) -> None
            sets the energy for this Snapshot

            C++ signature :
                void SetEnergy(RDKit::Snapshot {lvalue},double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, coordList: list, energy: float = 0.0) -> object: 
        """
        __init__( arg1: AtomPairsParameters, coordList: list, energy: float = 0.0) -> object
            Constructor;
            coordList: list of floats containing the coordinates for this Snapshot;
            energy:    the energy for this Snapshot.
            

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::list {lvalue} [,double=0.0])

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::Snapshot*)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, other: Snapshot) -> object: ...
    pass
class Trajectory(Boost.Python.instance):
    """
    A class which allows storing Snapshots from a trajectory
    """
    def AddConformersToMol(self, mol: Mol, fromCid: int = -1, toCid: int = -1) -> int: 
        """
        AddConformersToMol( self: Trajectory, mol: Mol, fromCid: int = -1, toCid: int = -1) -> int
            adds conformations from the Trajectory to mol
            fromCid is the first Snapshot that will be added as a Conformer; defaults to -1 (first available)
            toCid is the last Snapshot that will be added as a Conformer; defaults to -1 (all)
            

            C++ signature :
                unsigned int AddConformersToMol(RDKit::Trajectory {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1]])
        """
    def AddSnapshot(self, s: Snapshot) -> int: 
        """
        AddSnapshot( self: Trajectory, s: Snapshot) -> int
            appends Snapshot s to this Trajectory; returns the zero-based index position of the added snapshot
            

            C++ signature :
                unsigned int AddSnapshot(RDKit::Trajectory {lvalue},RDKit::Snapshot)
        """
    def Clear(self) -> None: 
        """
        Clear( self: Trajectory) -> None
            removes all Snapshots from the Trajectory
            

            C++ signature :
                void Clear(RDKit::Trajectory {lvalue})
        """
    def Dimension(self) -> int: 
        """
        Dimension( self: Trajectory) -> int
            returns the dimensionality of this Trajectory's coordinate tuples

            C++ signature :
                unsigned int Dimension(RDKit::Trajectory {lvalue})
        """
    def GetSnapshot(self, snapshotNum: int) -> Snapshot: 
        """
        GetSnapshot( self: Trajectory, snapshotNum: int) -> Snapshot
            returns the Snapshot snapshotNum, where the latter is the zero-based index of the retrieved Snapshot
            

            C++ signature :
                RDKit::Snapshot* GetSnapshot(RDKit::Trajectory*,unsigned int)
        """
    def InsertSnapshot(self, snapshotNum: int, s: Snapshot) -> int: 
        """
        InsertSnapshot( self: Trajectory, snapshotNum: int, s: Snapshot) -> int
            inserts Snapshot s into the Trajectory at the position snapshotNum, where the latter is the zero-based index of the Trajectory's Snapshot before which the Snapshot s will be inserted; returns the zero-based index position of the inserted snapshot
            

            C++ signature :
                unsigned int InsertSnapshot(RDKit::Trajectory {lvalue},unsigned int,RDKit::Snapshot)
        """
    def NumPoints(self) -> int: 
        """
        NumPoints( self: Trajectory) -> int
            returns the number of coordinate tuples associated to each Snapshot

            C++ signature :
                unsigned int NumPoints(RDKit::Trajectory {lvalue})
        """
    def RemoveSnapshot(self, snapshotNum: int) -> int: 
        """
        RemoveSnapshot( self: Trajectory, snapshotNum: int) -> int
            removes Snapshot snapshotNum from the Trajectory, where snapshotNum is the zero-based index of Snapshot to be removed
            

            C++ signature :
                unsigned int RemoveSnapshot(RDKit::Trajectory {lvalue},unsigned int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, dimension: int, numPoints: int, snapshotList: list = []) -> object: 
        """
        __init__( arg1: AtomPairsParameters, dimension: int, numPoints: int, snapshotList: list = []) -> object
            Constructor;
            dimension:    dimensionality of this Trajectory's coordinate tuples;
            numPoints:    number of coordinate tuples associated to each Snapshot;
            snapshotList: list of Snapshot objects used to initialize the Trajectory (optional; defaults to []).
            

            C++ signature :
                void* __init__(boost::python::api::object,unsigned int,unsigned int [,boost::python::list=[]])

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::Trajectory*)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, other: Trajectory) -> object: ...
    @staticmethod
    def __len__( arg1: Trajectory) -> int: 
        """
        __len__( arg1: Trajectory) -> int

            C++ signature :
                unsigned long __len__(RDKit::Trajectory {lvalue})
        """
    pass
def ReadAmberTrajectory( fName: str, traj: Trajectory) -> int:
    """
    ReadAmberTrajectory( fName: str, traj: Trajectory) -> int
        reads coordinates from an AMBER trajectory file into the Trajectory object; returns the number of Snapshot objects read in
        

        C++ signature :
            unsigned int ReadAmberTrajectory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::Trajectory {lvalue})
    """
def ReadGromosTrajectory( fName: str, traj: Trajectory) -> int:
    """
    ReadGromosTrajectory( fName: str, traj: Trajectory) -> int
        reads coordinates from a GROMOS trajectory file into the Trajectory object; returns the number of Snapshot objects read in
        

        C++ signature :
            unsigned int ReadGromosTrajectory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::Trajectory {lvalue})
    """
