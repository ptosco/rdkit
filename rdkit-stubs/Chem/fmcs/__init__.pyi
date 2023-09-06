""" Actual implementation of the FMCS algorithm
This code should be used by importing rdkit.Chem.MCS

"""
from __future__ import annotations
import rdkit.Chem.fmcs
import typing
from rdkit.Chem.fmcs.fmcs import Atom
from rdkit.Chem.fmcs.fmcs import AtomSmartsNoAromaticity
from rdkit.Chem.fmcs.fmcs import Bond
from rdkit.Chem.fmcs.fmcs import CachingTargetsMatcher
from rdkit.Chem.fmcs.fmcs import CangenNode
from collections import Counter
from rdkit.Chem.fmcs.fmcs import Default
from rdkit.Chem.fmcs.fmcs import DirectedEdge
from rdkit.Chem.fmcs.fmcs import FragmentedTypedMolecule
from rdkit.Chem.fmcs.fmcs import MCSResult
from rdkit.Chem.fmcs.fmcs import OutgoingEdge
from rdkit.Chem.fmcs.fmcs import SingleBestAtoms
from rdkit.Chem.fmcs.fmcs import SingleBestAtomsCompleteRingsOnly
from rdkit.Chem.fmcs.fmcs import SingleBestBonds
from rdkit.Chem.fmcs.fmcs import SingleBestBondsCompleteRingsOnly
from rdkit.Chem.fmcs.fmcs import Subgraph
from rdkit.Chem.fmcs.fmcs import Timer
from rdkit.Chem.fmcs.fmcs import TypedFragment
from rdkit.Chem.fmcs.fmcs import TypedMolecule
from rdkit.Chem.fmcs.fmcs import Uniquer
from rdkit.Chem.fmcs.fmcs import VerboseCachingTargetsMatcher
from rdkit.Chem.fmcs.fmcs import VerboseHeapOps
from itertools import chain
from itertools import combinations
from collections import defaultdict
from rdkit.Chem.fmcs.fmcs import starting_from
import copy
import itertools
import rdkit.Chem
import re
import sys
import time
import weakref

__all__ = [
    "Atom",
    "AtomSmartsNoAromaticity",
    "Bond",
    "CachingTargetsMatcher",
    "CangenNode",
    "Chem",
    "Counter",
    "Default",
    "DirectedEdge",
    "EnumerationMolecule",
    "FragmentedTypedMolecule",
    "MATCH",
    "MCSResult",
    "OutgoingEdge",
    "SingleBestAtoms",
    "SingleBestAtomsCompleteRingsOnly",
    "SingleBestBonds",
    "SingleBestBondsCompleteRingsOnly",
    "Subgraph",
    "Timer",
    "TypedFragment",
    "TypedMolecule",
    "Uniquer",
    "VerboseCachingTargetsMatcher",
    "VerboseHeapOps",
    "all_subgraph_extensions",
    "assign_isotopes_from_class_tag",
    "atom_typer_any",
    "atom_typer_elements",
    "atom_typer_isotopes",
    "atom_typers",
    "bond_typer_any",
    "bond_typer_bondtypes",
    "bond_typers",
    "canon",
    "chain",
    "check_completeRingsOnly",
    "combinations",
    "compare_shortcuts",
    "compute_mcs",
    "convert_input_to_typed_molecules",
    "copy",
    "default_atom_typer",
    "default_bond_typer",
    "defaultdict",
    "eleno",
    "enumerate_subgraphs",
    "find_duplicates",
    "find_extension_size",
    "find_extensions",
    "find_upper_fragment_size_limits",
    "fmcs",
    "fragmented_mol_to_enumeration_mols",
    "gen_primes",
    "generate_smarts",
    "get_canonical_bondtype_counts",
    "get_canonical_bondtypes",
    "get_closure_label",
    "get_counts",
    "get_initial_cangen_nodes",
    "get_isotopes",
    "get_selected_atom_classes",
    "get_specified_types",
    "get_typed_fragment",
    "get_typed_molecule",
    "heapify",
    "heappop",
    "heappush",
    "intersect_counts",
    "itertools",
    "main",
    "make_arbitrary_smarts",
    "make_canonical_smarts",
    "make_complete_sdf",
    "make_fragment_sdf",
    "make_fragment_smiles",
    "make_structure_format",
    "namedtuple",
    "nonempty_powerset",
    "parse_num_atoms",
    "parse_select",
    "parse_threshold",
    "parse_timeout",
    "powerset",
    "prune_maximize_atoms",
    "prune_maximize_bonds",
    "range_pat",
    "re",
    "remove_unknown_bondtypes",
    "rerank",
    "restore_isotopes",
    "save_atom_classes",
    "save_isotopes",
    "set_isotopes",
    "starting_from",
    "structure_format_functions",
    "subgraph_to_fragment",
    "sys",
    "tiebreaker",
    "time",
    "value_pat",
    "weakref"
]


atom_typers: dict # value = {'any': <function atom_typer_any>, 'elements': <function atom_typer_elements>, 'isotopes': <function atom_typer_isotopes>}
bond_typers: dict # value = {'any': <function bond_typer_any>, 'bondtypes': <function bond_typer_bondtypes>}
compare_shortcuts = {'topology': ('any', 'any'), 'elements': ('elements', 'any'), 'types': ('elements', 'bondtypes')}
eleno = 52
range_pat: re.Pattern # value = re.compile('(\\d+)-(\\d*)')
structure_format_functions: dict # value = {'fragment-smiles': <function make_fragment_smiles>, 'fragment-sdf': <function make_fragment_sdf>, 'complete-sdf': <function make_complete_sdf>}
value_pat: re.Pattern # value = re.compile('(\\d+)')
EnumerationMolecule = rdkit.Chem.fmcs.fmcs.Molecule
default_atom_typer = rdkit.Chem.fmcs.fmcs.atom_typer_elements
default_bond_typer = rdkit.Chem.fmcs.fmcs.bond_typer_bondtypes
tiebreaker = rdkit.Chem.fmcs.fmcs._Counter.<locals>.<lambda>
