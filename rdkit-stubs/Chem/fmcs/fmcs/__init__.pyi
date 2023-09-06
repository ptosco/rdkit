"""FMCS - Find Maximum Common Substructure

This software finds the maximum common substructure of a set of
structures and reports it as a SMARTS strings.

This implements what I think is a new algorithm for the MCS problem.
The core description is:

  best_substructure = None
  pick one structure as the query, and other as the targets
  for each substructure in the query graph:
    convert it to a SMARTS string based on the desired match properties
    if the SMARTS pattern exists in all of the targets:
       then this is a common substructure
       keep track of the maximum such common structure,

The SMARTS string depends on the desired match properties. For
example, if ring atoms are only allowed to match ring atoms then an
aliphatic ring carbon in the query is converted to the SMARTS "[C;R]",
and the double-bond ring bond converted to "=;@" while the respectice
chain-only version are "[C;!R]" and "=;!@".

The algorithm I outlined earlier will usually take a long time. There
are several ways to speed it up.

== Bond elimination ==

As the first step, remove bonds which obviously cannot be part of the
MCS.

This requires atom and bond type information, which I store as SMARTS
patterns. A bond can only be in the MCS if its canonical bond type is
present in all of the structures. A bond type is string made of the
SMARTS for one atom, the SMARTS for the bond, and the SMARTS for the
other atom. The canonical bond type is the lexographically smaller of
the two possible bond types for a bond.

The atom and bond SMARTS depend on the type comparison used.

The "ring-matches-ring-only" option adds an "@" or "!@" to the bond
SMARTS, so that the canonical bondtype for "C-C" becomes [#6]-@[#6] or
[#6]-!@[#6] if the bond is in a ring or not in a ring, and if atoms
are compared by element and bonds are compared by bondtype. (This
option does not add "R" or "!R" to the atom SMARTS because there
should be a single bond in the MCS of c1ccccc1O and CO.)

The result of all of this atom and bond typing is a "TypedMolecule"
for each input structure.

I then find which canonical bondtypes are present in all of the
structures. I convert each TypedMolecule into a
FragmentedTypedMolecule which has the same atom information but only
those bonds whose bondtypes are in all of the structures. This can
break a structure into multiple, disconnected fragments, hence the
name.

(BTW, I would like to use the fragmented molecules as the targets
because I think the SMARTS match would go faster, but the RDKit SMARTS
matcher doesn't like them. I think it's because the new molecule
hasn't been sanitized and the underlying data structure the ring
information doesn't exist. Instead, I use the input structures for the
SMARTS match.)

== Use the structure with the smallest largest fragment as the query ==
== and sort the targets by the smallest largest fragment             ==

I pick one of the FragmentedTypedMolecule instances as the source of
substructure enumeration. Which one?

My heuristic is to use the one with the smallest largest fragment.
Hopefully it produces the least number of subgraphs, but that's also
related to the number of rings, so a large linear graph will product
fewer subgraphs than a small fused ring system. I don't know how to
quantify that.

For each of the fragmented structures, I find the number of atoms in
the fragment with the most atoms, and I find the number of bonds in
the fragment with the most bonds. These might not be the same
fragment.

I sort the input structures by the number of bonds in the largest
fragment, with ties broken first on the number of atoms, and then on
the input order. The smallest such structure is the query structure,
and the remaining are the targets.

== Use a breadth-first search and a priority queue to    ==
== enumerate the fragment subgraphs                      ==

I extract each of the fragments from the FragmentedTypedMolecule into
a TypedFragment, which I use to make an EnumerationMolecule. An
enumeration molecule contains a pair of directed edges for each atom,
which simplifies the enumeration algorithm.

The enumeration algorithm is based around growing a seed. A seed
contains the current subgraph atoms and bonds as well as an exclusion
set of bonds which cannot be used for future grown. The initial seed
is the first bond in the fragment, which may potentially grow to use
the entire fragment. The second seed is the second bond in the
fragment, which is excluded from using the first bond in future
growth. The third seed starts from the third bond, which may not use
the first or second bonds during growth, and so on.


A seed can grow along bonds connected to an atom in the seed but which
aren't already in the seed and aren't in the set of excluded bonds for
the seed. If there are no such bonds then subgraph enumeration ends
for this fragment. Given N bonds there are 2**N-1 possible ways to
grow, which is just the powerset of the available bonds, excluding the
no-growth case.

This breadth-first growth takes into account all possibilties of using
the available N bonds so all of those bonds are added to the exclusion
set of the newly expanded subgraphs.

For performance reasons, the bonds used for growth are separated into
'internal' bonds, which connect two atoms already in the subgraph, and
'external' bonds, which lead outwards to an atom not already in the
subgraph.

Each seed growth can add from 0 to N new atoms and bonds. The goal is
to maximize the subgraph size so the seeds are stored in a priority
queue, ranked so the seed with the most bonds is processed first. This
turns the enumeration into something more like a depth-first search.


== Prune seeds which aren't found in all of the structures ==

At each stage of seed growth I check that the new seed exists in all
of the original structures. (Well, all except the one which I
enumerate over in the first place; by definition that one will match.)
If it doesn't match then there's no reason to include this seed or any
larger seeds made from it.

The check is easy; I turn the subgraph into its corresponding SMARTS
string and use RDKit's normal SMARTS matcher to test for a match.

There are three ways to generate a SMARTS string: 1) arbitrary, 2)
canonical, 3) hybrid.

I have not tested #1. During most of the development I assumed that
SMARTS matches across a few hundred structures would be slow, so that
the best solution is to generate a *canonical* SMARTS and cache the
match information.

Well, it turns out that my canonical SMARTS match code takes up most
of the FMCS run-time. If I drop the canonicalization step then the
code averages about 5-10% faster. This isn't the same as #1 - I still
do the initial atom assignment based on its neighborhood, which is
like a circular fingerprint of size 2 and *usually* gives a consistent
SMARTS pattern, which I can then cache.

However, there are times when the non-canonical SMARTS code is slower.
Obviously one is if there are a lot of structures, and another if is
there is a lot of symmetry. I'm still working on characterizing this.


== Maximize atoms? or bonds? ==

The above algorithm enumerates all subgraphs of the query and
identifies those subgraphs which are common to all input structures.

It's trivial then to keep track of the current "best" subgraph, which
can defined as having the subgraph with the most atoms, or the most
bonds. Both of those options are implemented.

It would not be hard to keep track of all other subgraphs which are
the same size.

== --complete-ring-only implementation ==

The "complete ring only" option is implemented by first enabling the
"ring-matches-ring-only" option, as otherwise it doesn't make sense.

Second, in order to be a "best" subgraph, all bonds in the subgraph
which are ring bonds in the original molecule must also be in a ring
in the subgraph. This is handled as a post-processing step.

(Note: some possible optimizations, like removing ring bonds from
structure fragments which are not in a ring, are not yet implemented.)


== Prune seeds which have no potential for growing large enough  ==

Given a seed, its set of edges available for growth, and the set of
excluded bonds, figure out the maximum possible growth for the seed.
If this maximum possible is less than the current best subgraph then
prune.

This requires a graph search, currently done in Python, which is a bit
expensive. To speed things up, I precompute some edge information.
That is, if I know that a given bond is a chain bond (not in a ring)
then I can calculate the maximum number of atoms and bonds for seed
growth along that bond, in either direction. However, precomputation
doesn't take into account the excluded bonds, so after a while the
predicted value is too high.

Again, I'm still working on characterizing this, and an implementation
in C++ would have different tradeoffs.
"""
from __future__ import annotations
import rdkit.Chem.fmcs.fmcs
import typing
from collections import Counter
from itertools import chain
from itertools import combinations
from collections import defaultdict
import _collections
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


class Atom(tuple):
    """
    Atom(real_atom, atom_smarts, bond_indices, is_in_ring)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('real_atom', 'atom_smarts', 'bond_indices', 'is_in_ring')
    _fields_defaults = {}
    atom_smarts: _collections._tuplegetter
    bond_indices: _collections._tuplegetter
    is_in_ring: _collections._tuplegetter
    real_atom: _collections._tuplegetter
    pass
class AtomSmartsNoAromaticity(dict):
    pass
class Bond(tuple):
    """
    Bond(real_bond, bond_smarts, canonical_bondtype, atom_indices, is_in_ring)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('real_bond', 'bond_smarts', 'canonical_bondtype', 'atom_indices', 'is_in_ring')
    _fields_defaults = {}
    atom_indices: _collections._tuplegetter
    bond_smarts: _collections._tuplegetter
    canonical_bondtype: _collections._tuplegetter
    is_in_ring: _collections._tuplegetter
    real_bond: _collections._tuplegetter
    pass
class CachingTargetsMatcher(dict):
    pass
class CangenNode():
    __slots__ = ['index', 'atom_smarts', 'value', 'neighbors', 'rank', 'outgoing_edges']
    atom_smarts: member_descriptor # value = <member 'atom_smarts' of 'CangenNode' objects>
    index: member_descriptor # value = <member 'index' of 'CangenNode' objects>
    neighbors: member_descriptor # value = <member 'neighbors' of 'CangenNode' objects>
    outgoing_edges: member_descriptor # value = <member 'outgoing_edges' of 'CangenNode' objects>
    rank: member_descriptor # value = <member 'rank' of 'CangenNode' objects>
    value: member_descriptor # value = <member 'value' of 'CangenNode' objects>
    pass
class Default():
    atomCompare = 'elements'
    bondCompare = 'bondtypes'
    completeRingsOnly = False
    matchValences = False
    maximize = 'bonds'
    ringMatchesRingOnly = False
    timeout = None
    timeoutString = 'none'
    pass
class DirectedEdge(tuple):
    """
    DirectedEdge(bond_index, end_atom_index)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('bond_index', 'end_atom_index')
    _fields_defaults = {}
    bond_index: _collections._tuplegetter
    end_atom_index: _collections._tuplegetter
    pass
class FragmentedTypedMolecule():
    pass
class MCSResult():
    pass
class OutgoingEdge(tuple):
    """
    OutgoingEdge(from_atom_index, bond_index, bond_smarts, other_node_idx, other_node)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('from_atom_index', 'bond_index', 'bond_smarts', 'other_node_idx', 'other_node')
    _fields_defaults = {}
    bond_index: _collections._tuplegetter
    bond_smarts: _collections._tuplegetter
    from_atom_index: _collections._tuplegetter
    other_node: _collections._tuplegetter
    other_node_idx: _collections._tuplegetter
    pass
class _SingleBest():
    pass
class SingleBestAtomsCompleteRingsOnly(_SingleBest):
    pass
class SingleBestBonds(_SingleBest):
    pass
class SingleBestBondsCompleteRingsOnly(_SingleBest):
    pass
class Subgraph(tuple):
    """
    Subgraph(atom_indices, bond_indices)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('atom_indices', 'bond_indices')
    _fields_defaults = {}
    atom_indices: _collections._tuplegetter
    bond_indices: _collections._tuplegetter
    pass
class Timer():
    pass
class TypedFragment():
    pass
class TypedMolecule():
    pass
class Uniquer(dict):
    pass
class VerboseCachingTargetsMatcher():
    pass
class VerboseHeapOps():
    pass
class SingleBestAtoms(_SingleBest):
    pass
class starting_from():
    pass
__version__ = '1.1'
__version_info = (1, 1, 0)
_atom_class_dict: weakref.WeakKeyDictionary # value = <WeakKeyDictionary>
_atom_smarts_no_aromaticity = {1: '#1', 5: '#5', 6: '#6', 7: '#7', 8: '#8', 15: '#15', 16: '#16', 33: '#33', 34: '#34', 52: '#52', 2: 'He'}
_available_closures = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
_get_symbol: method # value = <bound method GetElementSymbol of <rdkit.Chem.rdchem.PeriodicTable object>>
_isotope_dict: weakref.WeakKeyDictionary # value = <WeakKeyDictionary>
_maximize_options: dict # value = {('atoms', False): (<function prune_maximize_atoms>, <class 'rdkit.Chem.fmcs.fmcs.SingleBestAtoms'>), ('atoms', True): (<function prune_maximize_atoms>, <class 'rdkit.Chem.fmcs.fmcs.SingleBestAtomsCompleteRingsOnly'>), ('bonds', False): (<function prune_maximize_bonds>, <class 'rdkit.Chem.fmcs.fmcs.SingleBestBonds'>), ('bonds', True): (<function prune_maximize_bonds>, <class 'rdkit.Chem.fmcs.fmcs.SingleBestBondsCompleteRingsOnly'>)}
_prime_stream: generator # value = <generator object gen_primes>
_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927]
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
