""" making random numbers consistent so we get good regressions

"""
from __future__ import annotations
import rdkit.RDRandom
import typing
import random
import sys

__all__ = [
    "random",
    "randrange",
    "seed",
    "shuffle",
    "sys"
]


def random(*args, **kwargs) -> typing.Any:
    pass
randrange: method # value = <bound method Random.randrange of <random.Random object>>
seed: method # value = <bound method Random.seed of <random.Random object>>
shuffle: method # value = <bound method Random.shuffle of <random.Random object>>
