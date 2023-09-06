""" Automatic search for quantization bounds

This uses the expected informational gain to determine where quantization bounds should
lie.

**Notes**:

  - bounds are less than, so if the bounds are [1.,2.],
    [0.9,1.,1.1,2.,2.2] -> [0,1,1,2,2]

"""
from __future__ import annotations
import rdkit.ML.Data.Quantize
import typing
import numpy
import rdkit.ML.Data.cQuantize
import rdkit.ML.InfoTheory.entropy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "FindVarMultQuantBounds",
    "FindVarQuantBound",
    "cQuantize",
    "entropy",
    "feq",
    "hascQuantize",
    "numpy"
]


def _FindStartPoints( values: AtomPairsParameters, results: AtomPairsParameters, nData: int) -> list:
    """
    _FindStartPoints( values: AtomPairsParameters, results: AtomPairsParameters, nData: int) -> list
        TODO: provide docstring

        C++ signature :
            boost::python::list _FindStartPoints(boost::python::api::object,boost::python::api::object,int)
    """
def _RecurseOnBounds( vals: AtomPairsParameters, pyCuts: list, which: int, pyStarts: list, results: AtomPairsParameters, nPossibleRes: int) -> tuple:
    """
    _RecurseOnBounds( vals: AtomPairsParameters, pyCuts: list, which: int, pyStarts: list, results: AtomPairsParameters, nPossibleRes: int) -> tuple
        TODO: provide docstring

        C++ signature :
            boost::python::tuple _RecurseOnBounds(boost::python::api::object,boost::python::list,int,boost::python::list,boost::python::api::object,int)
    """
_float_tol = 1e-08
hascQuantize = 1
