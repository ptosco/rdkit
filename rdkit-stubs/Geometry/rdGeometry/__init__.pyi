"""Module containing geometry objects like points, grids, etc
"""
from __future__ import annotations
import rdkit.Geometry.rdGeometry
import typing
import Boost.Python

__all__ = [
    "ComputeDihedralAngle",
    "ComputeGridCentroid",
    "ComputeSignedDihedralAngle",
    "FindGridTerminalPoints",
    "Point2D",
    "Point3D",
    "PointND",
    "ProtrudeDistance",
    "TanimotoDistance",
    "TverskyIndex",
    "UniformGrid3D",
    "UniformGrid3D_",
    "WriteGridToFile"
]


class Point2D(Boost.Python.instance):
    """
    A class to represent a two-dimensional point
    """
    @staticmethod
    def AngleTo( arg1: Point2D, arg2: Point2D) -> float: 
        """
        AngleTo( arg1: Point2D, arg2: Point2D) -> float
            determines the angle between a vector to this point (between 0 and PI)

            C++ signature :
                double AngleTo(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def DirectionVector( arg1: Point2D, arg2: Point2D) -> Point2D: 
        """
        DirectionVector( arg1: Point2D, arg2: Point2D) -> Point2D
            return a normalized direction vector from this point to another

            C++ signature :
                RDGeom::Point2D DirectionVector(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def DotProduct( arg1: Point2D, arg2: Point2D) -> float: 
        """
        DotProduct( arg1: Point2D, arg2: Point2D) -> float
            Dot product with another point

            C++ signature :
                double DotProduct(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def Length( arg1: Point2D) -> float: 
        """
        Length( arg1: Point2D) -> float
            Length of the vector

            C++ signature :
                double Length(RDGeom::Point2D {lvalue})
        """
    @staticmethod
    def LengthSq( arg1: Point2D) -> float: 
        """
        LengthSq( arg1: Point2D) -> float
            Square of the length

            C++ signature :
                double LengthSq(RDGeom::Point2D {lvalue})
        """
    @staticmethod
    def Normalize( arg1: Point2D) -> None: 
        """
        Normalize( arg1: Point2D) -> None
            Normalize the vector (using L2 norm)

            C++ signature :
                void Normalize(RDGeom::Point2D {lvalue})
        """
    @staticmethod
    def SignedAngleTo( arg1: Point2D, arg2: Point2D) -> float: 
        """
        SignedAngleTo( arg1: Point2D, arg2: Point2D) -> float
            determines the signed angle between a vector to this point (between 0 and 2*PI)

            C++ signature :
                double SignedAngleTo(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def __add__( arg1: Point2D, arg2: Point2D) -> object: 
        """
        __add__( arg1: Point2D, arg2: Point2D) -> object

            C++ signature :
                _object* __add__(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def __getinitargs__( arg1: Point2D) -> tuple: 
        """
        __getinitargs__( arg1: Point2D) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDGeom::Point2D)
        """
    @staticmethod
    def __getitem__( arg1: Point2D, arg2: int) -> float: 
        """
        __getitem__( arg1: Point2D, arg2: int) -> float

            C++ signature :
                double __getitem__(RDGeom::Point2D,int)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    def __iadd__( arg1: object, arg2: Point2D) -> object: 
        """
        __iadd__( arg1: object, arg2: Point2D) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDGeom::Point2D&>,RDGeom::Point2D)
        """
    @staticmethod
    def __idiv__( arg1: Point2D, arg2: float) -> Point2D: 
        """
        __idiv__( arg1: Point2D, arg2: float) -> Point2D
            Scalar division

            C++ signature :
                RDGeom::Point2D {lvalue} __idiv__(RDGeom::Point2D {lvalue},double)
        """
    @staticmethod
    def __imul__( arg1: Point2D, arg2: float) -> Point2D: 
        """
        __imul__( arg1: Point2D, arg2: float) -> Point2D
            Scalar multiplication

            C++ signature :
                RDGeom::Point2D {lvalue} __imul__(RDGeom::Point2D {lvalue},double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None
            Default Constructor

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,double,double)

            C++ signature :
                void __init__(_object*,RDGeom::Point3D)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: float, arg3: float) -> None: ...
    @typing.overload
    def __init__(self, other: Point3D) -> None: ...
    @staticmethod
    def __isub__( arg1: object, arg2: Point2D) -> object: 
        """
        __isub__( arg1: object, arg2: Point2D) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDGeom::Point2D&>,RDGeom::Point2D)
        """
    @staticmethod
    def __len__( arg1: Point2D) -> int: 
        """
        __len__( arg1: Point2D) -> int

            C++ signature :
                unsigned int __len__(RDGeom::Point2D {lvalue})
        """
    @staticmethod
    def __mul__( arg1: Point2D, arg2: float) -> object: 
        """
        __mul__( arg1: Point2D, arg2: float) -> object

            C++ signature :
                _object* __mul__(RDGeom::Point2D {lvalue},double)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: Point2D, arg2: Point2D) -> object: 
        """
        __sub__( arg1: Point2D, arg2: Point2D) -> object

            C++ signature :
                _object* __sub__(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
    @staticmethod
    def __truediv__( arg1: Point2D, arg2: float) -> object: 
        """
        __truediv__( arg1: Point2D, arg2: float) -> object

            C++ signature :
                _object* __truediv__(RDGeom::Point2D {lvalue},double)
        """
    @property
    def x(self) -> None:
        """
        :type: None
        """
    @property
    def y(self) -> None:
        """
        :type: None
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 48
    __safe_for_unpickling__ = True
    pass
class Point3D(Boost.Python.instance):
    """
    A class to represent a three-dimensional point
    The x, y, and z coordinates can be read and written using either attributes
    (i.e. pt.x = 4) or indexing (i.e. pt[0] = 4).
    """
    @staticmethod
    def AngleTo( arg1: Point3D, arg2: Point3D) -> float: 
        """
        AngleTo( arg1: Point3D, arg2: Point3D) -> float
            determines the angle between a vector to this point (between 0 and PI)

            C++ signature :
                double AngleTo(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def CrossProduct( arg1: Point3D, arg2: Point3D) -> Point3D: 
        """
        CrossProduct( arg1: Point3D, arg2: Point3D) -> Point3D
            Get the cross product between two points

            C++ signature :
                RDGeom::Point3D CrossProduct(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def DirectionVector( arg1: Point3D, arg2: Point3D) -> Point3D: 
        """
        DirectionVector( arg1: Point3D, arg2: Point3D) -> Point3D
            return a normalized direction vector from this point to another

            C++ signature :
                RDGeom::Point3D DirectionVector(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def Distance( arg1: Point3D, arg2: Point3D) -> float: 
        """
        Distance( arg1: Point3D, arg2: Point3D) -> float
            Distance from this point to another point

            C++ signature :
                double Distance(RDGeom::Point3D,RDGeom::Point3D)
        """
    @staticmethod
    def DotProduct( arg1: Point3D, arg2: Point3D) -> float: 
        """
        DotProduct( arg1: Point3D, arg2: Point3D) -> float
            Dot product with another point

            C++ signature :
                double DotProduct(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def Length( arg1: Point3D) -> float: 
        """
        Length( arg1: Point3D) -> float
            Length of the vector

            C++ signature :
                double Length(RDGeom::Point3D {lvalue})
        """
    @staticmethod
    def LengthSq( arg1: Point3D) -> float: 
        """
        LengthSq( arg1: Point3D) -> float
            Square of the length

            C++ signature :
                double LengthSq(RDGeom::Point3D {lvalue})
        """
    @staticmethod
    def Normalize( arg1: Point3D) -> None: 
        """
        Normalize( arg1: Point3D) -> None
            Normalize the vector (using L2 norm)

            C++ signature :
                void Normalize(RDGeom::Point3D {lvalue})
        """
    @staticmethod
    def SignedAngleTo( arg1: Point3D, arg2: Point3D) -> float: 
        """
        SignedAngleTo( arg1: Point3D, arg2: Point3D) -> float
            determines the signed angle between a vector to this point (between 0 and 2*PI)

            C++ signature :
                double SignedAngleTo(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def __add__( arg1: Point3D, arg2: Point3D) -> object: 
        """
        __add__( arg1: Point3D, arg2: Point3D) -> object

            C++ signature :
                _object* __add__(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def __getinitargs__( arg1: Point3D) -> tuple: 
        """
        __getinitargs__( arg1: Point3D) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDGeom::Point3D)
        """
    @staticmethod
    def __getitem__( arg1: Point3D, arg2: int) -> float: 
        """
        __getitem__( arg1: Point3D, arg2: int) -> float

            C++ signature :
                double __getitem__(RDGeom::Point3D,int)
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
    def __iadd__( arg1: Point3D, arg2: Point3D) -> Point3D: 
        """
        __iadd__( arg1: Point3D, arg2: Point3D) -> Point3D
            Addition to another point

            C++ signature :
                RDGeom::Point3D {lvalue} __iadd__(RDGeom::Point3D {lvalue},RDGeom::Point3D)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDGeom::Point3D&>,RDGeom::Point3D)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: Point3D) -> object: ...
    @staticmethod
    def __idiv__( arg1: Point3D, arg2: float) -> Point3D: 
        """
        __idiv__( arg1: Point3D, arg2: float) -> Point3D
            Scalar division

            C++ signature :
                RDGeom::Point3D {lvalue} __idiv__(RDGeom::Point3D {lvalue},double)
        """
    @staticmethod
    def __imul__( arg1: Point3D, arg2: float) -> Point3D: 
        """
        __imul__( arg1: Point3D, arg2: float) -> Point3D
            Scalar multiplication

            C++ signature :
                RDGeom::Point3D {lvalue} __imul__(RDGeom::Point3D {lvalue},double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None
            Default Constructor

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,double,double,double)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: float, arg3: float, arg4: float) -> None: ...
    @staticmethod
    @typing.overload
    def __isub__( arg1: Point3D, arg2: Point3D) -> Point3D: 
        """
        __isub__( arg1: Point3D, arg2: Point3D) -> Point3D
            Vector difference

            C++ signature :
                RDGeom::Point3D {lvalue} __isub__(RDGeom::Point3D {lvalue},RDGeom::Point3D)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDGeom::Point3D&>,RDGeom::Point3D)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: Point3D) -> object: ...
    @staticmethod
    def __len__( arg1: Point3D) -> int: 
        """
        __len__( arg1: Point3D) -> int

            C++ signature :
                unsigned int __len__(RDGeom::Point3D {lvalue})
        """
    @staticmethod
    def __mul__( arg1: Point3D, arg2: float) -> object: 
        """
        __mul__( arg1: Point3D, arg2: float) -> object

            C++ signature :
                _object* __mul__(RDGeom::Point3D {lvalue},double)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: Point3D, arg2: Point3D) -> object: 
        """
        __sub__( arg1: Point3D, arg2: Point3D) -> object

            C++ signature :
                _object* __sub__(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def __truediv__( arg1: Point3D, arg2: float) -> object: 
        """
        __truediv__( arg1: Point3D, arg2: float) -> object

            C++ signature :
                _object* __truediv__(RDGeom::Point3D {lvalue},double)
        """
    @property
    def x(self) -> None:
        """
        :type: None
        """
    @property
    def y(self) -> None:
        """
        :type: None
        """
    @property
    def z(self) -> None:
        """
        :type: None
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 56
    __safe_for_unpickling__ = True
    pass
class PointND(Boost.Python.instance):
    """
    A class to represent an N-dimensional point
    """
    @staticmethod
    def AngleTo( arg1: PointND, arg2: PointND) -> float: 
        """
        AngleTo( arg1: PointND, arg2: PointND) -> float
            determines the angle between a vector to this point (between 0 and PI)

            C++ signature :
                double AngleTo(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
    @staticmethod
    def DirectionVector( arg1: PointND, arg2: PointND) -> PointND: 
        """
        DirectionVector( arg1: PointND, arg2: PointND) -> PointND
            return a normalized direction vector from this point to another

            C++ signature :
                RDGeom::PointND DirectionVector(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
    @staticmethod
    def Distance( arg1: Point3D, arg2: Point3D) -> float: 
        """
        Distance( arg1: Point3D, arg2: Point3D) -> float
            Distance from this point to another point

            C++ signature :
                double Distance(RDGeom::Point3D,RDGeom::Point3D)
        """
    @staticmethod
    def DotProduct( arg1: PointND, arg2: PointND) -> float: 
        """
        DotProduct( arg1: PointND, arg2: PointND) -> float
            Dot product with another point

            C++ signature :
                double DotProduct(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
    @staticmethod
    def Length( arg1: PointND) -> float: 
        """
        Length( arg1: PointND) -> float
            Length of the vector

            C++ signature :
                double Length(RDGeom::PointND {lvalue})
        """
    @staticmethod
    def LengthSq( arg1: PointND) -> float: 
        """
        LengthSq( arg1: PointND) -> float
            Square of the length

            C++ signature :
                double LengthSq(RDGeom::PointND {lvalue})
        """
    @staticmethod
    def Normalize( arg1: PointND) -> None: 
        """
        Normalize( arg1: PointND) -> None
            Normalize the vector (using L2 norm)

            C++ signature :
                void Normalize(RDGeom::PointND {lvalue})
        """
    @staticmethod
    def __add__( arg1: PointND, arg2: PointND) -> object: 
        """
        __add__( arg1: PointND, arg2: PointND) -> object

            C++ signature :
                _object* __add__(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
    @staticmethod
    def __getinitargs__( arg1: PointND) -> tuple: 
        """
        __getinitargs__( arg1: PointND) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDGeom::PointND)
        """
    @staticmethod
    def __getitem__( arg1: PointND, arg2: int) -> float: 
        """
        __getitem__( arg1: PointND, arg2: int) -> float

            C++ signature :
                double __getitem__(RDGeom::PointND,int)
        """
    @staticmethod
    def __getstate__( arg1: PointND) -> tuple: 
        """
        __getstate__( arg1: PointND) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(RDGeom::PointND)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: PointND, arg2: PointND) -> PointND: 
        """
        __iadd__( arg1: PointND, arg2: PointND) -> PointND
            Addition to another point

            C++ signature :
                RDGeom::PointND {lvalue} __iadd__(RDGeom::PointND {lvalue},RDGeom::PointND)

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDGeom::PointND&>,RDGeom::PointND)
        """
    @staticmethod
    @typing.overload
    def __iadd__( arg1: object, arg2: PointND) -> object: ...
    @staticmethod
    def __idiv__( arg1: PointND, arg2: float) -> PointND: 
        """
        __idiv__( arg1: PointND, arg2: float) -> PointND
            Scalar division

            C++ signature :
                RDGeom::PointND {lvalue} __idiv__(RDGeom::PointND {lvalue},double)
        """
    @staticmethod
    def __imul__( arg1: PointND, arg2: float) -> PointND: 
        """
        __imul__( arg1: PointND, arg2: float) -> PointND
            Scalar multiplication

            C++ signature :
                RDGeom::PointND {lvalue} __imul__(RDGeom::PointND {lvalue},double)
        """
    @staticmethod
    def __init__( arg1: object, arg2: int) -> None: 
        """
        __init__( arg1: object, arg2: int) -> None

            C++ signature :
                void __init__(_object*,unsigned int)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: PointND, arg2: PointND) -> PointND: 
        """
        __isub__( arg1: PointND, arg2: PointND) -> PointND
            Vector difference

            C++ signature :
                RDGeom::PointND {lvalue} __isub__(RDGeom::PointND {lvalue},RDGeom::PointND)

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDGeom::PointND&>,RDGeom::PointND)
        """
    @staticmethod
    @typing.overload
    def __isub__( arg1: object, arg2: PointND) -> object: ...
    @staticmethod
    def __len__( arg1: PointND) -> int: 
        """
        __len__( arg1: PointND) -> int

            C++ signature :
                unsigned int __len__(RDGeom::PointND {lvalue})
        """
    @staticmethod
    def __mul__( arg1: PointND, arg2: float) -> object: 
        """
        __mul__( arg1: PointND, arg2: float) -> object

            C++ signature :
                _object* __mul__(RDGeom::PointND {lvalue},double)
        """
    @staticmethod
    def __setitem__( arg1: PointND, arg2: int, arg3: float) -> float: 
        """
        __setitem__( arg1: PointND, arg2: int, arg3: float) -> float

            C++ signature :
                double __setitem__(RDGeom::PointND {lvalue},int,double)
        """
    @staticmethod
    def __setstate__( arg1: PointND, arg2: tuple) -> None: 
        """
        __setstate__( arg1: PointND, arg2: tuple) -> None

            C++ signature :
                void __setstate__(RDGeom::PointND {lvalue},boost::python::tuple)
        """
    @staticmethod
    def __sub__( arg1: PointND, arg2: PointND) -> object: 
        """
        __sub__( arg1: PointND, arg2: PointND) -> object

            C++ signature :
                _object* __sub__(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
    @staticmethod
    def __truediv__( arg1: PointND, arg2: float) -> object: 
        """
        __truediv__( arg1: PointND, arg2: float) -> object

            C++ signature :
                _object* __truediv__(RDGeom::PointND {lvalue},double)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 48
    __safe_for_unpickling__ = True
    pass
class UniformGrid3D_(Boost.Python.instance):
    """
    Class to represent a uniform three-dimensional
        cubic grid. Each grid point can store a poisitive integer value. For the sake
        of efficiency these value can either be binary, fit in 2, 4, 8 or 16 bits
    """
    @staticmethod
    def CompareParams( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> bool: 
        """
        CompareParams( arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> bool
            Compare the parameters between two grid object

            C++ signature :
                bool CompareParams(RDGeom::UniformGrid3D {lvalue},RDGeom::UniformGrid3D)
        """
    @staticmethod
    def GetGridIndex( arg1: UniformGrid3D_, arg2: int, arg3: int, arg4: int) -> int: 
        """
        GetGridIndex( arg1: UniformGrid3D_, arg2: int, arg3: int, arg4: int) -> int
            Get the index to the grid point with the three integer indices provided

            C++ signature :
                int GetGridIndex(RDGeom::UniformGrid3D {lvalue},unsigned int,unsigned int,unsigned int)
        """
    @staticmethod
    def GetGridIndices( arg1: UniformGrid3D_, arg2: int) -> tuple: 
        """
        GetGridIndices( arg1: UniformGrid3D_, arg2: int) -> tuple
            Returns the integer indices of the grid index provided.

            C++ signature :
                boost::python::tuple GetGridIndices(RDGeom::UniformGrid3D,unsigned int)
        """
    @staticmethod
    def GetGridPointIndex( arg1: UniformGrid3D_, arg2: Point3D) -> int: 
        """
        GetGridPointIndex( arg1: UniformGrid3D_, arg2: Point3D) -> int
            Get the index to the grid point closest to the specified point

            C++ signature :
                int GetGridPointIndex(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D)
        """
    @staticmethod
    def GetGridPointLoc( arg1: UniformGrid3D_, arg2: int) -> Point3D: 
        """
        GetGridPointLoc( arg1: UniformGrid3D_, arg2: int) -> Point3D
            Get the location of the specified grid point

            C++ signature :
                RDGeom::Point3D GetGridPointLoc(RDGeom::UniformGrid3D {lvalue},unsigned int)
        """
    @staticmethod
    def GetNumX( arg1: UniformGrid3D_) -> int: 
        """
        GetNumX( arg1: UniformGrid3D_) -> int
            Get the number of grid points along x-axis

            C++ signature :
                unsigned int GetNumX(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetNumY( arg1: UniformGrid3D_) -> int: 
        """
        GetNumY( arg1: UniformGrid3D_) -> int
            Get the number of grid points along y-axis

            C++ signature :
                unsigned int GetNumY(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetNumZ( arg1: UniformGrid3D_) -> int: 
        """
        GetNumZ( arg1: UniformGrid3D_) -> int
            Get the number of grid points along z-axis

            C++ signature :
                unsigned int GetNumZ(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetOccupancyVect( arg1: UniformGrid3D_) -> DiscreteValueVect: 
        """
        GetOccupancyVect( arg1: UniformGrid3D_) -> DiscreteValueVect
            Get the occupancy vector for the grid

            C++ signature :
                RDKit::DiscreteValueVect const* GetOccupancyVect(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetOffset( arg1: UniformGrid3D_) -> Point3D: 
        """
        GetOffset( arg1: UniformGrid3D_) -> Point3D
            Get the location of the center of the grid

            C++ signature :
                RDGeom::Point3D GetOffset(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetSize( arg1: UniformGrid3D_) -> int: 
        """
        GetSize( arg1: UniformGrid3D_) -> int
            Get the size of the grid (number of grid points)

            C++ signature :
                unsigned int GetSize(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetSpacing( arg1: UniformGrid3D_) -> float: 
        """
        GetSpacing( arg1: UniformGrid3D_) -> float
            Get the grid spacing

            C++ signature :
                double GetSpacing(RDGeom::UniformGrid3D {lvalue})
        """
    @staticmethod
    def GetVal( arg1: UniformGrid3D_, arg2: int) -> int: 
        """
        GetVal( arg1: UniformGrid3D_, arg2: int) -> int
            Get the value at the specified grid point

            C++ signature :
                int GetVal(RDGeom::UniformGrid3D,unsigned int)
        """
    @staticmethod
    def GetValPoint( arg1: UniformGrid3D_, arg2: Point3D) -> int: 
        """
        GetValPoint( arg1: UniformGrid3D_, arg2: Point3D) -> int
            Get the value at the closest grid point

            C++ signature :
                int GetValPoint(RDGeom::UniformGrid3D,RDGeom::Point3D)
        """
    def SetSphereOccupancy(self, center: Point3D, radius: float, stepSize: float, maxLayers: int = -1, ignoreOutOfBound: bool = True) -> None: 
        """
        SetSphereOccupancy( self: UniformGrid3D_, center: Point3D, radius: float, stepSize: float, maxLayers: int = -1, ignoreOutOfBound: bool = True) -> None
            Set the occupancy on the grid for a sphere or specified radius
             and multiple layers around this sphere, with decreasing values of 
            occupancy
            

            C++ signature :
                void SetSphereOccupancy(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D,double,double [,int=-1 [,bool=True]])
        """
    @staticmethod
    def SetVal( arg1: UniformGrid3D_, arg2: int, arg3: int) -> None: 
        """
        SetVal( arg1: UniformGrid3D_, arg2: int, arg3: int) -> None
            Set the value at the specified grid point

            C++ signature :
                void SetVal(RDGeom::UniformGrid3D {lvalue},unsigned int,unsigned int)
        """
    @staticmethod
    def SetValPoint( arg1: UniformGrid3D_, arg2: Point3D, arg3: int) -> None: 
        """
        SetValPoint( arg1: UniformGrid3D_, arg2: Point3D, arg3: int) -> None
            Set the value at grid point closest to the specified point

            C++ signature :
                void SetValPoint(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D,unsigned int)
        """
    @staticmethod
    def __getinitargs__( arg1: UniformGrid3D_) -> tuple: 
        """
        __getinitargs__( arg1: UniformGrid3D_) -> tuple

            C++ signature :
                boost::python::tuple __getinitargs__(RDGeom::UniformGrid3D)
        """
    @staticmethod
    def __getstate__( arg1: AtomPairsParameters) -> tuple: 
        """
        __getstate__( arg1: AtomPairsParameters) -> tuple

            C++ signature :
                boost::python::tuple __getstate__(boost::python::api::object)
        """
    @staticmethod
    def __iadd__( arg1: object, arg2: UniformGrid3D_) -> object: 
        """
        __iadd__( arg1: object, arg2: UniformGrid3D_) -> object

            C++ signature :
                _object* __iadd__(boost::python::back_reference<RDGeom::UniformGrid3D&>,RDGeom::UniformGrid3D)
        """
    @staticmethod
    def __iand__( arg1: object, arg2: UniformGrid3D_) -> object: 
        """
        __iand__( arg1: object, arg2: UniformGrid3D_) -> object

            C++ signature :
                _object* __iand__(boost::python::back_reference<RDGeom::UniformGrid3D&>,RDGeom::UniformGrid3D)
        """
    @staticmethod
    def __init__( arg1: object, arg2: str) -> None: 
        """
        __init__( arg1: object, arg2: str) -> None
            pickle constructor

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __ior__( arg1: object, arg2: UniformGrid3D_) -> object: 
        """
        __ior__( arg1: object, arg2: UniformGrid3D_) -> object

            C++ signature :
                _object* __ior__(boost::python::back_reference<RDGeom::UniformGrid3D&>,RDGeom::UniformGrid3D)
        """
    @staticmethod
    def __isub__( arg1: object, arg2: UniformGrid3D_) -> object: 
        """
        __isub__( arg1: object, arg2: UniformGrid3D_) -> object

            C++ signature :
                _object* __isub__(boost::python::back_reference<RDGeom::UniformGrid3D&>,RDGeom::UniformGrid3D)
        """
    @staticmethod
    def __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None: 
        """
        __setstate__( arg1: AtomPairsParameters, arg2: tuple) -> None

            C++ signature :
                void __setstate__(boost::python::api::object,boost::python::tuple)
        """
    __getstate_manages_dict__ = True
    __instance_size__ = 96
    __safe_for_unpickling__ = True
    pass
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
