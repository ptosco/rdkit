""" functions to match a bunch of fragment descriptors from a file

No user-servicable parts inside.  ;-)

"""
from __future__ import annotations
import rdkit.Chem.Fragments
import typing
import os
import rdkit.Chem
import rdkit.RDConfig

__all__ = [
    "Chem",
    "RDConfig",
    "defaultPatternFileName",
    "fn",
    "fns",
    "fr_Al_COO",
    "fr_Al_OH",
    "fr_Al_OH_noTert",
    "fr_ArN",
    "fr_Ar_COO",
    "fr_Ar_N",
    "fr_Ar_NH",
    "fr_Ar_OH",
    "fr_COO",
    "fr_COO2",
    "fr_C_O",
    "fr_C_O_noCOO",
    "fr_C_S",
    "fr_HOCCN",
    "fr_Imine",
    "fr_NH0",
    "fr_NH1",
    "fr_NH2",
    "fr_N_O",
    "fr_Ndealkylation1",
    "fr_Ndealkylation2",
    "fr_Nhpyrrole",
    "fr_SH",
    "fr_aldehyde",
    "fr_alkyl_carbamate",
    "fr_alkyl_halide",
    "fr_allylic_oxid",
    "fr_amide",
    "fr_amidine",
    "fr_aniline",
    "fr_aryl_methyl",
    "fr_azide",
    "fr_azo",
    "fr_barbitur",
    "fr_benzene",
    "fr_benzodiazepine",
    "fr_bicyclic",
    "fr_diazo",
    "fr_dihydropyridine",
    "fr_epoxide",
    "fr_ester",
    "fr_ether",
    "fr_furan",
    "fr_guanido",
    "fr_halogen",
    "fr_hdrzine",
    "fr_hdrzone",
    "fr_imidazole",
    "fr_imide",
    "fr_isocyan",
    "fr_isothiocyan",
    "fr_ketone",
    "fr_ketone_Topliss",
    "fr_lactam",
    "fr_lactone",
    "fr_methoxy",
    "fr_morpholine",
    "fr_nitrile",
    "fr_nitro",
    "fr_nitro_arom",
    "fr_nitro_arom_nonortho",
    "fr_nitroso",
    "fr_oxazole",
    "fr_oxime",
    "fr_para_hydroxylation",
    "fr_phenol",
    "fr_phenol_noOrthoHbond",
    "fr_phos_acid",
    "fr_phos_ester",
    "fr_piperdine",
    "fr_piperzine",
    "fr_priamide",
    "fr_prisulfonamd",
    "fr_pyridine",
    "fr_quatN",
    "fr_sulfide",
    "fr_sulfonamd",
    "fr_sulfone",
    "fr_term_acetylene",
    "fr_tetrazole",
    "fr_thiazole",
    "fr_thiocyan",
    "fr_thiophene",
    "fr_unbrch_alkane",
    "fr_urea",
    "name",
    "os"
]


defaultPatternFileName = '/scratch/toscopa1/src/rdkit/Data/FragmentDescriptors.csv'
fn = None
fns: list # value = [('fr_C_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_O_noCOO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_methoxy', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxime', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ether', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aldehyde', <function _LoadPatterns.<locals>.<lambda>>), ('fr_quatN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH0', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_N', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_NH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aniline', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Imine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitrile', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitroso', <function _LoadPatterns.<locals>.<lambda>>), ('fr_N_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_diazo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_priamide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amidine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_guanido', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Nhpyrrole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isothiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_halogen', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_halide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_SH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_S', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_prisulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_barbitur', <function _LoadPatterns.<locals>.<lambda>>), ('fr_urea', <function _LoadPatterns.<locals>.<lambda>>), ('fr_term_acetylene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imidazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_furan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiophene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_pyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperdine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_morpholine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactam', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_tetrazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_epoxide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_unbrch_alkane', <function _LoadPatterns.<locals>.<lambda>>), ('fr_bicyclic', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_acid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom_nonortho', <function _LoadPatterns.<locals>.<lambda>>), ('fr_dihydropyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol_noOrthoHbond', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH_noTert', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzodiazepine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_para_hydroxylation', <function _LoadPatterns.<locals>.<lambda>>), ('fr_allylic_oxid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aryl_methyl', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_carbamate', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone_Topliss', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ArN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_HOCCN', <function _LoadPatterns.<locals>.<lambda>>)]
name = 'fr_HOCCN'
fr_Al_COO = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Al_OH = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Al_OH_noTert = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_ArN = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ar_COO = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ar_N = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ar_NH = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ar_OH = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_COO = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_COO2 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_C_O = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_C_O_noCOO = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_C_S = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_HOCCN = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Imine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_NH0 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_NH1 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_NH2 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_N_O = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ndealkylation1 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Ndealkylation2 = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_Nhpyrrole = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_SH = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_aldehyde = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_alkyl_carbamate = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_alkyl_halide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_allylic_oxid = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_amide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_amidine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_aniline = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_aryl_methyl = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_azide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_azo = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_barbitur = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_benzene = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_benzodiazepine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_bicyclic = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_diazo = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_dihydropyridine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_epoxide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_ester = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_ether = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_furan = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_guanido = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_halogen = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_hdrzine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_hdrzone = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_imidazole = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_imide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_isocyan = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_isothiocyan = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_ketone = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_ketone_Topliss = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_lactam = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_lactone = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_methoxy = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_morpholine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_nitrile = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_nitro = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_nitro_arom = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_nitro_arom_nonortho = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_nitroso = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_oxazole = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_oxime = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_para_hydroxylation = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_phenol = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_phenol_noOrthoHbond = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_phos_acid = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_phos_ester = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_piperdine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_piperzine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_priamide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_prisulfonamd = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_pyridine = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_quatN = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_sulfide = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_sulfonamd = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_sulfone = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_term_acetylene = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_tetrazole = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_thiazole = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_thiocyan = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_thiophene = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_unbrch_alkane = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
fr_urea = rdkit.Chem.Fragments._LoadPatterns.<locals>.<lambda>
