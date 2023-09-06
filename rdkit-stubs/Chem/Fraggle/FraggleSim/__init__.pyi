"""
Fragmentation algorithm
-----------------------

identify acyclic bonds
enumerate all single cuts
make sure you chop off more that 1 atom
keeps bits which are >60% query mol
enumerate all double cuts
keeps bits with 1 attachment point (i.e throw middle bit away)
need to be >60% query mol

identify exocyclic bonds
enumerate all single "ring" cuts
Check if it results in more that one component
keep correct bit if >40% query mol

enumerate successful "rings" cuts with an acyclic cut
Check if it results in more that one component
keep correct if >60% query mol

"""
from __future__ import annotations
import rdkit.Chem.Fraggle.FraggleSim
import typing
from itertools import combinations
import rdkit.Chem
import rdkit.Chem.rdchem
import rdkit.Chem.rdqueries
import rdkit.DataStructs
import sys

__all__ = [
    "ACYC_SMARTS",
    "CYC_SMARTS",
    "Chem",
    "DataStructs",
    "FTYPE_ACYCLIC",
    "FTYPE_CYCLIC",
    "FTYPE_CYCLIC_ACYCLIC",
    "GetFraggleSimilarity",
    "atomContrib",
    "cSma1",
    "cSma2",
    "combinations",
    "compute_fraggle_similarity_for_subs",
    "delete_bonds",
    "dummyAtomQuery",
    "generate_fraggle_fragmentation",
    "isValidRingCut",
    "modified_query_fps",
    "rdkitFpParams",
    "rdqueries",
    "select_fragments",
    "sys"
]


ACYC_SMARTS: rdkit.Chem.rdchem.Mol
CYC_SMARTS: rdkit.Chem.rdchem.Mol
FTYPE_ACYCLIC = 'acyclic'
FTYPE_CYCLIC = 'cyclic'
FTYPE_CYCLIC_ACYCLIC = 'cyclic_and_acyclic'
cSma1: rdkit.Chem.rdchem.Mol
cSma2: rdkit.Chem.rdchem.Mol
dummyAtomQuery: rdkit.Chem.rdchem.QueryAtom
modified_query_fps = {}
rdkitFpParams = {'maxPath': 5, 'fpSize': 1024, 'nBitsPerHash': 2}
