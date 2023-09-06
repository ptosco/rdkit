""" Calculation of Lipinski parameters for molecules

"""
from __future__ import annotations
import rdkit.Chem.Lipinski
import typing
import rdkit.Chem
import rdkit.Chem.rdMolDescriptors
import rdkit.Chem.rdchem

__all__ = [
    "Chem",
    "FractionCSP3",
    "HAcceptorSmarts",
    "HDonorSmarts",
    "HeavyAtomCount",
    "HeteroatomSmarts",
    "NHOHCount",
    "NHOHSmarts",
    "NOCount",
    "NOCountSmarts",
    "NumAliphaticCarbocycles",
    "NumAliphaticHeterocycles",
    "NumAliphaticRings",
    "NumAromaticCarbocycles",
    "NumAromaticHeterocycles",
    "NumAromaticRings",
    "NumHAcceptors",
    "NumHDonors",
    "NumHeteroatoms",
    "NumRotatableBonds",
    "NumSaturatedCarbocycles",
    "NumSaturatedHeterocycles",
    "NumSaturatedRings",
    "RingCount",
    "RotatableBondSmarts",
    "nm",
    "rdMolDescriptors",
    "txt"
]


HAcceptorSmarts: rdkit.Chem.rdchem.Mol
HDonorSmarts: rdkit.Chem.rdchem.Mol
HeteroatomSmarts: rdkit.Chem.rdchem.Mol
NHOHSmarts: rdkit.Chem.rdchem.Mol
NOCountSmarts: rdkit.Chem.rdchem.Mol
RotatableBondSmarts: rdkit.Chem.rdchem.Mol
_bulkConvert = ('CalcFractionCSP3', 'CalcNumAromaticRings', 'CalcNumSaturatedRings', 'CalcNumAromaticHeterocycles', 'CalcNumAromaticCarbocycles', 'CalcNumSaturatedHeterocycles', 'CalcNumSaturatedCarbocycles', 'CalcNumAliphaticRings', 'CalcNumAliphaticHeterocycles', 'CalcNumAliphaticCarbocycles')
nm = 'NumAliphaticCarbocycles'
txt = 'CalcNumAliphaticCarbocycles'
FractionCSP3 = rdkit.Chem.Lipinski.<lambda>
NHOHCount = rdkit.Chem.Lipinski.<lambda>
NOCount = rdkit.Chem.Lipinski.<lambda>
NumAliphaticCarbocycles = rdkit.Chem.Lipinski.<lambda>
NumAliphaticHeterocycles = rdkit.Chem.Lipinski.<lambda>
NumAliphaticRings = rdkit.Chem.Lipinski.<lambda>
NumAromaticCarbocycles = rdkit.Chem.Lipinski.<lambda>
NumAromaticHeterocycles = rdkit.Chem.Lipinski.<lambda>
NumAromaticRings = rdkit.Chem.Lipinski.<lambda>
NumHAcceptors = rdkit.Chem.Lipinski.<lambda>
NumHDonors = rdkit.Chem.Lipinski.<lambda>
NumHeteroatoms = rdkit.Chem.Lipinski.<lambda>
NumRotatableBonds = rdkit.Chem.Lipinski.<lambda>
NumSaturatedCarbocycles = rdkit.Chem.Lipinski.<lambda>
NumSaturatedHeterocycles = rdkit.Chem.Lipinski.<lambda>
NumSaturatedRings = rdkit.Chem.Lipinski.<lambda>
RingCount = rdkit.Chem.Lipinski.<lambda>
_HAcceptors = rdkit.Chem.Lipinski.<lambda>
_HDonors = rdkit.Chem.Lipinski.<lambda>
_Heteroatoms = rdkit.Chem.Lipinski.<lambda>
_RotatableBonds = rdkit.Chem.Lipinski.<lambda>
_fn = rdkit.Chem.Lipinski.<lambda>
