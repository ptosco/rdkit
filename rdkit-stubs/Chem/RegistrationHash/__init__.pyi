"""
Generate a unique hash code for a molecule based on chemistry. If two
molecules are chemically "the same", they should have the same hash.

Using molhash adds value beyond using SMILES because it:

* Ignores SMILES features that are not chemically meaningful
(e.g. atom map numbers)
* Canonicalizes enhanced stereochemistry groups. For example
`C[C@H](O)CC |&1:1|` and `C[C@@H](O)CC |&1:1|` have the same
molhash
* Canonicalizes S group data (for example, polymer data)

There are two hash schemes, the default, and one in which
tautomers are considered equivalent.

"""
from __future__ import annotations
import rdkit.Chem.RegistrationHash
import typing
import enum
import hashlib
import json
import logging
import rdkit.Chem
import rdkit.Chem.rdMolHash
import re
import types

__all__ = [
    "ATOM_PROP_MAP_NUMBER",
    "Chem",
    "DEFAULT_CXFLAG",
    "EMPTY_MOL_TAUTOMER_HASH",
    "ENHANCED_STEREO_GROUP_REGEX",
    "ENHANCED_STEREO_GROUP_WEIGHTS",
    "EnhancedStereoUpdateMode",
    "GetMolHash",
    "GetMolLayers",
    "GetNoStereoLayers",
    "GetStereoTautomerHash",
    "HashLayer",
    "HashScheme",
    "Iterable",
    "Optional",
    "enum",
    "hashlib",
    "json",
    "logger",
    "logging",
    "rdMolHash",
    "re"
]


class EnhancedStereoUpdateMode(enum.Enum):
    """
    An enumeration.
    """
    ADD_WEIGHTS: rdkit.Chem.RegistrationHash.EnhancedStereoUpdateMode # value = <EnhancedStereoUpdateMode.ADD_WEIGHTS: 1>
    REMOVE_WEIGHTS: rdkit.Chem.RegistrationHash.EnhancedStereoUpdateMode # value = <EnhancedStereoUpdateMode.REMOVE_WEIGHTS: 2>
    __members__: mappingproxy # value = mappingproxy({'ADD_WEIGHTS': <EnhancedStereoUpdateMode.ADD_WEIGHTS: 1>, 'REMOVE_WEIGHTS': <EnhancedStereoUpdateMode.REMOVE_WEIGHTS: 2>})
    name: types.DynamicClassAttribute
    value: types.DynamicClassAttribute
    pass
class HashLayer(enum.Enum):
    """
    :cvar CANONICAL_SMILES: RDKit canonical SMILES (excluding enhanced stereo)
    :cvar ESCAPE: arbitrary other information to be incorporated
    :cvar FORMULA: a simple molecular formula for the molecule
    :cvar NO_STEREO_SMILES: RDKit canonical SMILES with all stereo removed
    :cvar SGROUP_DATA: canonicalization of all SGroups data present
    :cvar TAUTOMER_HASH: SMILES-like representation for a generic tautomer form
    :cvar NO_STEREO_TAUTOMER_HASH: the above tautomer hash lacking all stereo
    """
    CANONICAL_SMILES: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.CANONICAL_SMILES: 1>
    ESCAPE: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.ESCAPE: 2>
    FORMULA: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.FORMULA: 3>
    NO_STEREO_SMILES: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.NO_STEREO_SMILES: 4>
    NO_STEREO_TAUTOMER_HASH: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>
    SGROUP_DATA: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.SGROUP_DATA: 6>
    TAUTOMER_HASH: rdkit.Chem.RegistrationHash.HashLayer # value = <HashLayer.TAUTOMER_HASH: 7>
    __members__: mappingproxy # value = mappingproxy({'CANONICAL_SMILES': <HashLayer.CANONICAL_SMILES: 1>, 'ESCAPE': <HashLayer.ESCAPE: 2>, 'FORMULA': <HashLayer.FORMULA: 3>, 'NO_STEREO_SMILES': <HashLayer.NO_STEREO_SMILES: 4>, 'NO_STEREO_TAUTOMER_HASH': <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, 'SGROUP_DATA': <HashLayer.SGROUP_DATA: 6>, 'TAUTOMER_HASH': <HashLayer.TAUTOMER_HASH: 7>})
    name: types.DynamicClassAttribute
    value: types.DynamicClassAttribute
    pass
class HashScheme(enum.Enum):
    """
    Which hash layers to use to when deduplicating molecules

    Typically the "ALL_LAYERS" scheme is used, but some users may want
    the "TAUTOMER_INSENSITIVE_LAYERS" scheme.

    :cvar ALL_LAYERS: most strict hash scheme utilizing all layers
    :cvar STEREO_INSENSITIVE_LAYERS: excludes stereo sensitive layers
    :cvar TAUTOMER_INSENSITIVE_LAYERS: excludes tautomer sensitive layers
    """
    ALL_LAYERS: rdkit.Chem.RegistrationHash.HashScheme # value = <HashScheme.ALL_LAYERS: (<HashLayer.CANONICAL_SMILES: 1>, <HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_SMILES: 4>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>, <HashLayer.TAUTOMER_HASH: 7>)>
    STEREO_INSENSITIVE_LAYERS: rdkit.Chem.RegistrationHash.HashScheme # value = <HashScheme.STEREO_INSENSITIVE_LAYERS: (<HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_SMILES: 4>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>)>
    TAUTOMER_INSENSITIVE_LAYERS: rdkit.Chem.RegistrationHash.HashScheme # value = <HashScheme.TAUTOMER_INSENSITIVE_LAYERS: (<HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>, <HashLayer.TAUTOMER_HASH: 7>)>
    __members__: mappingproxy # value = mappingproxy({'ALL_LAYERS': <HashScheme.ALL_LAYERS: (<HashLayer.CANONICAL_SMILES: 1>, <HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_SMILES: 4>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>, <HashLayer.TAUTOMER_HASH: 7>)>, 'STEREO_INSENSITIVE_LAYERS': <HashScheme.STEREO_INSENSITIVE_LAYERS: (<HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_SMILES: 4>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>)>, 'TAUTOMER_INSENSITIVE_LAYERS': <HashScheme.TAUTOMER_INSENSITIVE_LAYERS: (<HashLayer.ESCAPE: 2>, <HashLayer.FORMULA: 3>, <HashLayer.NO_STEREO_TAUTOMER_HASH: 5>, <HashLayer.SGROUP_DATA: 6>, <HashLayer.TAUTOMER_HASH: 7>)>})
    name: types.DynamicClassAttribute
    value: types.DynamicClassAttribute
    pass
ATOM_PROP_MAP_NUMBER = 'molAtomMapNumber'
DEFAULT_CXFLAG = 65
EMPTY_MOL_TAUTOMER_HASH = '_0_0'
ENHANCED_STEREO_GROUP_REGEX: re.Pattern # value = re.compile('((?:a|[&o]\\d+):\\d+(?:,\\d+)*)')
ENHANCED_STEREO_GROUP_WEIGHTS = {rdkit.Chem.rdchem.StereoGroupType.STEREO_AND: 1000, rdkit.Chem.rdchem.StereoGroupType.STEREO_OR: 2000, rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE: 3000}
Iterable: typing._GenericAlias # value = typing.Iterable
Optional: typing._SpecialForm # value = typing.Optional
logger: logging.Logger # value = <Logger rdkit.Chem.RegistrationHash (WARNING)>
