""" A module for Geometry stuff	
 
"""
from __future__ import annotations
import rdkit.Geometry
import typing
from rdkit.Geometry.rdGeometry import Point2D
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Geometry.rdGeometry import PointND
from rdkit.Geometry.rdGeometry import UniformGrid3D_
import rdkit.DataStructs

__all__ = [
    "ComputeDihedralAngle",
    "ComputeGridCentroid",
    "ComputeSignedDihedralAngle",
    "DataStructs",
    "FindGridTerminalPoints",
    "Point2D",
    "Point3D",
    "PointND",
    "ProtrudeDistance",
    "TanimotoDistance",
    "TverskyIndex",
    "UniformGrid3D",
    "UniformGrid3D_",
    "WriteGridToFile",
    "rdGeometry"
]


def ComputeDihedralAngle( arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D) -> float:
    """
    ComputeDihedralAngle( arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D) -> float
        calculates the dihedral angle determined by four Point3D objects

        C++ signature :
            double ComputeDihedralAngle(RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D)
    """
def ComputeGridCentroid( arg1: UniformGrid3D_, arg2: Point3D, arg3: float) -> tuple:
    """
    ComputeGridCentroid( arg1: UniformGrid3D_, arg2: Point3D, arg3: float) -> tuple
        Compute the grid point at the center of sphere around a Point3D

        C++ signature :
            boost::python::tuple ComputeGridCentroid(RDGeom::UniformGrid3D,RDGeom::Point3D,double)
    """
def ComputeSignedDihedralAngle( arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D) -> float:
    """
    ComputeSignedDihedralAngle( arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D) -> float
        calculates the signed dihedral angle determined by four Point3D objects

        C++ signature :
            double ComputeSignedDihedralAngle(RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D)
    """
def FindGridTerminalPoints( arg1: UniformGrid3D_, arg2: float, arg3: float) -> tuple:
    """
    FindGridTerminalPoints( arg1: UniformGrid3D_, arg2: float, arg3: float) -> tuple
        Find a grid's terminal points (defined in the subshape algorithm).

        C++ signature :
            boost::python::tuple FindGridTerminalPoints(RDGeom::UniformGrid3D,double,double)
    """
def ProtrudeDistance( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float:
    """
    ProtrudeDistance( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float
        Compute the protrude distance between two grid objects

        C++ signature :
            double ProtrudeDistance(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D)
    """
def TanimotoDistance( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float:
    """
    TanimotoDistance( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float
        Compute the tanimoto distance between two grid objects

        C++ signature :
            double TanimotoDistance(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D)
    """
def TverskyIndex( arg1: UniformGrid3D_, arg2: UniformGrid3D_, arg3: float, arg4: float) -> float:
    """
    TverskyIndex( arg1: UniformGrid3D_, arg2: UniformGrid3D_, arg3: float, arg4: float) -> float
        Compute the tversky index between two grid objects

        C++ signature :
            double TverskyIndex(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D,double,double)
    """
def UniformGrid3D( dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, valType: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, offSet: Point3D = None) -> UniformGrid3D_:
    """
    UniformGrid3D( dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, valType: DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, offSet: Point3D = None) -> UniformGrid3D_
        Faking the constructor

        C++ signature :
            RDGeom::UniformGrid3D* UniformGrid3D(double,double,double [,double=0.5 [,RDKit::DiscreteValueVect::DiscreteValueType=rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE [,RDGeom::Point3D const*=None]]])
    """
def WriteGridToFile( arg1: UniformGrid3D_, arg2: str) -> None:
    """
    WriteGridToFile( arg1: UniformGrid3D_, arg2: str) -> None
        Write the grid to a grid file

        C++ signature :
            void WriteGridToFile(RDGeom::UniformGrid3D,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
