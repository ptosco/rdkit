from __future__ import annotations
import rdkit.Chem.MolKey.MolKey
import typing
import _collections
import base64
import hashlib
import logging
import os
import rdkit.Avalon.pyAvalonTools
import rdkit.Chem
import rdkit.Chem.MolKey.InchiInfo
import rdkit.RDConfig
import re
import tempfile
import uuid

__all__ = [
    "BAD_SET",
    "BadMoleculeException",
    "CHIRAL_POS",
    "CheckCTAB",
    "Chem",
    "ERROR_DICT",
    "ErrorBitsToText",
    "GET_STEREO_RE",
    "GetInchiForCTAB",
    "GetKeyForCTAB",
    "INCHI_COMPUTATION_ERROR",
    "INCHI_READWRITE_ERROR",
    "InchiInfo",
    "InchiResult",
    "MOL_KEY_VERSION",
    "MolIdentifierException",
    "MolKeyResult",
    "NULL_MOL",
    "NULL_SMILES_RE",
    "PATTERN_NULL_MOL",
    "RDConfig",
    "RDKIT_CONVERSION_ERROR",
    "T_NULL_MOL",
    "base64",
    "hashlib",
    "initStruchk",
    "logging",
    "namedtuple",
    "os",
    "pyAvalonTools",
    "re",
    "stereo_code_dict",
    "tempfile",
    "uuid"
]


class BadMoleculeException(Exception, BaseException):
    pass
class InchiResult(tuple):
    """
    InchiResult(error, inchi, fixed_ctab)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('error', 'inchi', 'fixed_ctab')
    _fields_defaults = {}
    error: _collections._tuplegetter
    fixed_ctab: _collections._tuplegetter
    inchi: _collections._tuplegetter
    pass
class MolIdentifierException(Exception, BaseException):
    pass
class MolKeyResult(tuple):
    """
    MolKeyResult(mol_key, error, inchi, fixed_ctab, stereo_code, stereo_comment)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('mol_key', 'error', 'inchi', 'fixed_ctab', 'stereo_code', 'stereo_comment')
    _fields_defaults = {}
    error: _collections._tuplegetter
    fixed_ctab: _collections._tuplegetter
    inchi: _collections._tuplegetter
    mol_key: _collections._tuplegetter
    stereo_code: _collections._tuplegetter
    stereo_comment: _collections._tuplegetter
    pass
BAD_SET = 986019
CHIRAL_POS = 12
ERROR_DICT = {'BAD_MOLECULE': 1, 'ALIAS_CONVERSION_FAILED': 2, 'TRANSFORMED': 4, 'FRAGMENTS_FOUND': 8, 'EITHER_WARNING': 16, 'STEREO_ERROR': 32, 'DUBIOUS_STEREO_REMOVED': 64, 'ATOM_CLASH': 128, 'ATOM_CHECK_FAILED': 256, 'SIZE_CHECK_FAILED': 512, 'RECHARGED': 1024, 'STEREO_FORCED_BAD': 2048, 'STEREO_TRANSFORMED': 4096, 'TEMPLATE_TRANSFORMED': 8192, 'INCHI_COMPUTATION_ERROR': 65536, 'RDKIT_CONVERSION_ERROR': 131072, 'INCHI_READWRITE_ERROR': 262144, 'NULL_MOL': 524288}
GET_STEREO_RE: re.Pattern # value = re.compile('^InChI=1S(.*?)/(t.*?)/m\\d/s1(.*$)')
INCHI_COMPUTATION_ERROR = 65536
INCHI_READWRITE_ERROR = 262144
MOL_KEY_VERSION = '1'
NULL_MOL = 524288
NULL_SMILES_RE: re.Pattern # value = re.compile('^\\s*$|^\\s*NO_STRUCTURE\\s*$', re.IGNORECASE)
PATTERN_NULL_MOL = '^([\\s0]+[1-9]+[\\s]+V[\\w]*)'
RDKIT_CONVERSION_ERROR = 131072
T_NULL_MOL = (524288, '')
__initCalled = False
stereo_code_dict = {'DEFAULT': 0, 'S_ACHIR': 1, 'S_ABS': 2, 'S_REL': 3, 'S_PART': 4, 'S_UNKN': 5, 'S_ABS_ACHIR': 6, 'R_ONE': 11, 'R_REL': 12, 'R_OTHER': 13, 'MX_ENANT': 21, 'MX_DIAST': 22, 'MX_SP2': 31, 'MX_DIAST_ABS': 32, 'MX_DIAST_REL': 33, 'OTHER': 100, 'UNDEFINED': 200}
