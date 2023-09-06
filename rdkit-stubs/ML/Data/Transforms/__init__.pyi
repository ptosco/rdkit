from __future__ import annotations
import rdkit.ML.Data.Transforms
import typing

__all__ = [
    "GetAvailTransforms"
]


_availTransforms: list # value = [('Center', <function _CenterTForm>, 'translates so that mean(x)=0'), ('Normalize', <function _NormalizeTForm>, 'scales so that dot(x,x)=1'), ('Standardize', <function _StandardTForm>, 'scales so that dev(x)=0')]
