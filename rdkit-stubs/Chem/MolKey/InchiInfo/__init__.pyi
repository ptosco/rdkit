from __future__ import annotations
import rdkit.Chem.MolKey.InchiInfo
import typing
import logging
import rdkit.Chem
import rdkit.Chem.inchi
import re

__all__ = [
    "Chem",
    "InchiInfo",
    "UPD_APP",
    "all_stereo_re",
    "console",
    "defined_stereo_re",
    "fixed_h_re",
    "h_layer_re",
    "inchi",
    "isotope_re",
    "logging",
    "mobile_h_atoms_re",
    "mobile_h_group_re",
    "re",
    "reconnected_re",
    "stereo_all_re",
    "stereo_re",
    "undef_stereo_re",
    "version_re"
]


class InchiInfo():
    pass
UPD_APP: logging.Logger # value = <Logger inchiinfo.application (INFO)>
all_stereo_re: re.Pattern # value = re.compile('(\\d+)[?+-]')
console: logging.StreamHandler # value = <StreamHandler <stderr> (NOTSET)>
defined_stereo_re: re.Pattern # value = re.compile('(\\d+)[+-]')
fixed_h_re: re.Pattern # value = re.compile('(.*?)/f(.*)')
h_layer_re: re.Pattern # value = re.compile('.*/h(.*)/?')
isotope_re: re.Pattern # value = re.compile('(.*?)/i(.*)')
mobile_h_atoms_re: re.Pattern # value = re.compile(',(\\d+)')
mobile_h_group_re: re.Pattern # value = re.compile('(\\(H.+?\\))')
reconnected_re: re.Pattern # value = re.compile('(.*?)/r(.*)')
stereo_all_re: re.Pattern # value = re.compile('.*/t([^/]+)')
stereo_re: re.Pattern # value = re.compile('.*/t(.*?)/.*')
undef_stereo_re: re.Pattern # value = re.compile('(\\d+)\\?')
version_re: re.Pattern # value = re.compile('(.*?)/(.*)')
