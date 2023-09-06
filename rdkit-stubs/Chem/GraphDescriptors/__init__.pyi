""" Calculation of topological/topochemical descriptors.



"""
from __future__ import annotations
import rdkit.Chem.GraphDescriptors
import typing
import math
import numpy
import rdkit.Chem
import rdkit.Chem.Graphs
import rdkit.Chem.rdMolDescriptors
import rdkit.Chem.rdchem
import rdkit.ML.InfoTheory.entropy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AvgIpc",
    "BalabanJ",
    "BertzCT",
    "Chem",
    "Chi0",
    "Chi0n",
    "Chi0v",
    "Chi1",
    "Chi1n",
    "Chi1v",
    "Chi2n",
    "Chi2v",
    "Chi3n",
    "Chi3v",
    "Chi4n",
    "Chi4v",
    "ChiNn_",
    "ChiNv_",
    "Graphs",
    "HallKierAlpha",
    "Ipc",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "entropy",
    "hallKierAlphas",
    "math",
    "numpy",
    "ptable",
    "rdMolDescriptors",
    "rdchem"
]


_log2val = 0.6931471805599453
hallKierAlphas = {'Br': [None, None, 0.48], 'C': [-0.22, -0.13, 0.0], 'Cl': [None, None, 0.29], 'F': [None, None, -0.07], 'H': [0.0, 0.0, 0.0], 'I': [None, None, 0.73], 'N': [-0.29, -0.2, -0.04], 'O': [None, -0.2, -0.04], 'P': [None, 0.3, 0.43], 'S': [None, 0.22, 0.35]}
ptable: rdkit.Chem.rdchem.PeriodicTable
Chi0n = rdkit.Chem.GraphDescriptors.<lambda>
Chi0v = rdkit.Chem.GraphDescriptors.<lambda>
Chi1n = rdkit.Chem.GraphDescriptors.<lambda>
Chi1v = rdkit.Chem.GraphDescriptors.<lambda>
Chi2n = rdkit.Chem.GraphDescriptors.<lambda>
Chi2v = rdkit.Chem.GraphDescriptors.<lambda>
Chi3n = rdkit.Chem.GraphDescriptors.<lambda>
Chi3v = rdkit.Chem.GraphDescriptors.<lambda>
Chi4n = rdkit.Chem.GraphDescriptors.<lambda>
Chi4v = rdkit.Chem.GraphDescriptors.<lambda>
ChiNn_ = rdkit.Chem.GraphDescriptors.<lambda>
ChiNv_ = rdkit.Chem.GraphDescriptors.<lambda>
HallKierAlpha = rdkit.Chem.GraphDescriptors.<lambda>
Kappa1 = rdkit.Chem.GraphDescriptors.<lambda>
Kappa2 = rdkit.Chem.GraphDescriptors.<lambda>
Kappa3 = rdkit.Chem.GraphDescriptors.<lambda>
