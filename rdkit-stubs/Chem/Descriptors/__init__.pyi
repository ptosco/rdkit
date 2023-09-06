from __future__ import annotations
import rdkit.Chem.Descriptors
import typing
import Boost.Python
import collections.abc
import rdkit.Chem
import rdkit.Chem.ChemUtils.DescriptorUtilities
import rdkit.Chem.rdMolDescriptors
import rdkit.Chem.rdPartialCharges

__all__ = [
    "AUTOCORR2D_1",
    "AUTOCORR2D_10",
    "AUTOCORR2D_100",
    "AUTOCORR2D_101",
    "AUTOCORR2D_102",
    "AUTOCORR2D_103",
    "AUTOCORR2D_104",
    "AUTOCORR2D_105",
    "AUTOCORR2D_106",
    "AUTOCORR2D_107",
    "AUTOCORR2D_108",
    "AUTOCORR2D_109",
    "AUTOCORR2D_11",
    "AUTOCORR2D_110",
    "AUTOCORR2D_111",
    "AUTOCORR2D_112",
    "AUTOCORR2D_113",
    "AUTOCORR2D_114",
    "AUTOCORR2D_115",
    "AUTOCORR2D_116",
    "AUTOCORR2D_117",
    "AUTOCORR2D_118",
    "AUTOCORR2D_119",
    "AUTOCORR2D_12",
    "AUTOCORR2D_120",
    "AUTOCORR2D_121",
    "AUTOCORR2D_122",
    "AUTOCORR2D_123",
    "AUTOCORR2D_124",
    "AUTOCORR2D_125",
    "AUTOCORR2D_126",
    "AUTOCORR2D_127",
    "AUTOCORR2D_128",
    "AUTOCORR2D_129",
    "AUTOCORR2D_13",
    "AUTOCORR2D_130",
    "AUTOCORR2D_131",
    "AUTOCORR2D_132",
    "AUTOCORR2D_133",
    "AUTOCORR2D_134",
    "AUTOCORR2D_135",
    "AUTOCORR2D_136",
    "AUTOCORR2D_137",
    "AUTOCORR2D_138",
    "AUTOCORR2D_139",
    "AUTOCORR2D_14",
    "AUTOCORR2D_140",
    "AUTOCORR2D_141",
    "AUTOCORR2D_142",
    "AUTOCORR2D_143",
    "AUTOCORR2D_144",
    "AUTOCORR2D_145",
    "AUTOCORR2D_146",
    "AUTOCORR2D_147",
    "AUTOCORR2D_148",
    "AUTOCORR2D_149",
    "AUTOCORR2D_15",
    "AUTOCORR2D_150",
    "AUTOCORR2D_151",
    "AUTOCORR2D_152",
    "AUTOCORR2D_153",
    "AUTOCORR2D_154",
    "AUTOCORR2D_155",
    "AUTOCORR2D_156",
    "AUTOCORR2D_157",
    "AUTOCORR2D_158",
    "AUTOCORR2D_159",
    "AUTOCORR2D_16",
    "AUTOCORR2D_160",
    "AUTOCORR2D_161",
    "AUTOCORR2D_162",
    "AUTOCORR2D_163",
    "AUTOCORR2D_164",
    "AUTOCORR2D_165",
    "AUTOCORR2D_166",
    "AUTOCORR2D_167",
    "AUTOCORR2D_168",
    "AUTOCORR2D_169",
    "AUTOCORR2D_17",
    "AUTOCORR2D_170",
    "AUTOCORR2D_171",
    "AUTOCORR2D_172",
    "AUTOCORR2D_173",
    "AUTOCORR2D_174",
    "AUTOCORR2D_175",
    "AUTOCORR2D_176",
    "AUTOCORR2D_177",
    "AUTOCORR2D_178",
    "AUTOCORR2D_179",
    "AUTOCORR2D_18",
    "AUTOCORR2D_180",
    "AUTOCORR2D_181",
    "AUTOCORR2D_182",
    "AUTOCORR2D_183",
    "AUTOCORR2D_184",
    "AUTOCORR2D_185",
    "AUTOCORR2D_186",
    "AUTOCORR2D_187",
    "AUTOCORR2D_188",
    "AUTOCORR2D_189",
    "AUTOCORR2D_19",
    "AUTOCORR2D_190",
    "AUTOCORR2D_191",
    "AUTOCORR2D_192",
    "AUTOCORR2D_2",
    "AUTOCORR2D_20",
    "AUTOCORR2D_21",
    "AUTOCORR2D_22",
    "AUTOCORR2D_23",
    "AUTOCORR2D_24",
    "AUTOCORR2D_25",
    "AUTOCORR2D_26",
    "AUTOCORR2D_27",
    "AUTOCORR2D_28",
    "AUTOCORR2D_29",
    "AUTOCORR2D_3",
    "AUTOCORR2D_30",
    "AUTOCORR2D_31",
    "AUTOCORR2D_32",
    "AUTOCORR2D_33",
    "AUTOCORR2D_34",
    "AUTOCORR2D_35",
    "AUTOCORR2D_36",
    "AUTOCORR2D_37",
    "AUTOCORR2D_38",
    "AUTOCORR2D_39",
    "AUTOCORR2D_4",
    "AUTOCORR2D_40",
    "AUTOCORR2D_41",
    "AUTOCORR2D_42",
    "AUTOCORR2D_43",
    "AUTOCORR2D_44",
    "AUTOCORR2D_45",
    "AUTOCORR2D_46",
    "AUTOCORR2D_47",
    "AUTOCORR2D_48",
    "AUTOCORR2D_49",
    "AUTOCORR2D_5",
    "AUTOCORR2D_50",
    "AUTOCORR2D_51",
    "AUTOCORR2D_52",
    "AUTOCORR2D_53",
    "AUTOCORR2D_54",
    "AUTOCORR2D_55",
    "AUTOCORR2D_56",
    "AUTOCORR2D_57",
    "AUTOCORR2D_58",
    "AUTOCORR2D_59",
    "AUTOCORR2D_6",
    "AUTOCORR2D_60",
    "AUTOCORR2D_61",
    "AUTOCORR2D_62",
    "AUTOCORR2D_63",
    "AUTOCORR2D_64",
    "AUTOCORR2D_65",
    "AUTOCORR2D_66",
    "AUTOCORR2D_67",
    "AUTOCORR2D_68",
    "AUTOCORR2D_69",
    "AUTOCORR2D_7",
    "AUTOCORR2D_70",
    "AUTOCORR2D_71",
    "AUTOCORR2D_72",
    "AUTOCORR2D_73",
    "AUTOCORR2D_74",
    "AUTOCORR2D_75",
    "AUTOCORR2D_76",
    "AUTOCORR2D_77",
    "AUTOCORR2D_78",
    "AUTOCORR2D_79",
    "AUTOCORR2D_8",
    "AUTOCORR2D_80",
    "AUTOCORR2D_81",
    "AUTOCORR2D_82",
    "AUTOCORR2D_83",
    "AUTOCORR2D_84",
    "AUTOCORR2D_85",
    "AUTOCORR2D_86",
    "AUTOCORR2D_87",
    "AUTOCORR2D_88",
    "AUTOCORR2D_89",
    "AUTOCORR2D_9",
    "AUTOCORR2D_90",
    "AUTOCORR2D_91",
    "AUTOCORR2D_92",
    "AUTOCORR2D_93",
    "AUTOCORR2D_94",
    "AUTOCORR2D_95",
    "AUTOCORR2D_96",
    "AUTOCORR2D_97",
    "AUTOCORR2D_98",
    "AUTOCORR2D_99",
    "AvgIpc",
    "BCUT2D_CHGHI",
    "BCUT2D_CHGLO",
    "BCUT2D_LOGPHI",
    "BCUT2D_LOGPLOW",
    "BCUT2D_MRHI",
    "BCUT2D_MRLOW",
    "BCUT2D_MWHI",
    "BCUT2D_MWLOW",
    "BalabanJ",
    "BertzCT",
    "CalcMolDescriptors",
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
    "EState_VSA1",
    "EState_VSA10",
    "EState_VSA11",
    "EState_VSA2",
    "EState_VSA3",
    "EState_VSA4",
    "EState_VSA5",
    "EState_VSA6",
    "EState_VSA7",
    "EState_VSA8",
    "EState_VSA9",
    "ExactMolWt",
    "FpDensityMorgan1",
    "FpDensityMorgan2",
    "FpDensityMorgan3",
    "FractionCSP3",
    "HallKierAlpha",
    "HeavyAtomCount",
    "HeavyAtomMolWt",
    "Ipc",
    "Kappa1",
    "Kappa2",
    "Kappa3",
    "LabuteASA",
    "MaxAbsEStateIndex",
    "MaxAbsPartialCharge",
    "MaxEStateIndex",
    "MaxPartialCharge",
    "MinAbsEStateIndex",
    "MinAbsPartialCharge",
    "MinEStateIndex",
    "MinPartialCharge",
    "MolLogP",
    "MolMR",
    "MolWt",
    "NHOHCount",
    "NOCount",
    "NumAliphaticCarbocycles",
    "NumAliphaticHeterocycles",
    "NumAliphaticRings",
    "NumAromaticCarbocycles",
    "NumAromaticHeterocycles",
    "NumAromaticRings",
    "NumHAcceptors",
    "NumHDonors",
    "NumHeteroatoms",
    "NumRadicalElectrons",
    "NumRotatableBonds",
    "NumSaturatedCarbocycles",
    "NumSaturatedHeterocycles",
    "NumSaturatedRings",
    "NumValenceElectrons",
    "PEOE_VSA1",
    "PEOE_VSA10",
    "PEOE_VSA11",
    "PEOE_VSA12",
    "PEOE_VSA13",
    "PEOE_VSA14",
    "PEOE_VSA2",
    "PEOE_VSA3",
    "PEOE_VSA4",
    "PEOE_VSA5",
    "PEOE_VSA6",
    "PEOE_VSA7",
    "PEOE_VSA8",
    "PEOE_VSA9",
    "PropertyFunctor",
    "RingCount",
    "SMR_VSA1",
    "SMR_VSA10",
    "SMR_VSA2",
    "SMR_VSA3",
    "SMR_VSA4",
    "SMR_VSA5",
    "SMR_VSA6",
    "SMR_VSA7",
    "SMR_VSA8",
    "SMR_VSA9",
    "SlogP_VSA1",
    "SlogP_VSA10",
    "SlogP_VSA11",
    "SlogP_VSA12",
    "SlogP_VSA2",
    "SlogP_VSA3",
    "SlogP_VSA4",
    "SlogP_VSA5",
    "SlogP_VSA6",
    "SlogP_VSA7",
    "SlogP_VSA8",
    "SlogP_VSA9",
    "TPSA",
    "VSA_EState1",
    "VSA_EState10",
    "VSA_EState2",
    "VSA_EState3",
    "VSA_EState4",
    "VSA_EState5",
    "VSA_EState6",
    "VSA_EState7",
    "VSA_EState8",
    "VSA_EState9",
    "abc",
    "autocorr",
    "descList",
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
    "names",
    "qed",
    "rdMolDescriptors",
    "rdPartialCharges",
    "setupAUTOCorrDescriptors"
]


class PropertyFunctor(rdkit.Chem.rdMolDescriptors.PythonPropertyFunctor, Boost.Python.instance):
    """
    Creates a python based property function that can be added to the
        global property list.  To use, subclass this class and override the
        __call__ method.  Then create an instance and add it to the
        registry.  The __call__ method should return a numeric value.

        Example:

          class NumAtoms(Descriptors.PropertyFunctor):
            def __init__(self):
              Descriptors.PropertyFunctor.__init__(self, "NumAtoms", "1.0.0")
            def __call__(self, mol):
              return mol.GetNumAtoms()

          numAtoms = NumAtoms()
          rdMolDescriptors.Properties.RegisterProperty(numAtoms)
        
    """
    pass
_descList: list # value = [('MaxAbsEStateIndex', <function MaxAbsEStateIndex>), ('MaxEStateIndex', <function MaxEStateIndex>), ('MinAbsEStateIndex', <function MinAbsEStateIndex>), ('MinEStateIndex', <function MinEStateIndex>), ('qed', <function qed>), ('MolWt', <function <lambda>>), ('HeavyAtomMolWt', <function HeavyAtomMolWt>), ('ExactMolWt', <function <lambda>>), ('NumValenceElectrons', <function NumValenceElectrons>), ('NumRadicalElectrons', <function NumRadicalElectrons>), ('MaxPartialCharge', <function MaxPartialCharge>), ('MinPartialCharge', <function MinPartialCharge>), ('MaxAbsPartialCharge', <function MaxAbsPartialCharge>), ('MinAbsPartialCharge', <function MinAbsPartialCharge>), ('FpDensityMorgan1', <function FpDensityMorgan1>), ('FpDensityMorgan2', <function FpDensityMorgan2>), ('FpDensityMorgan3', <function FpDensityMorgan3>), ('BCUT2D_MWHI', <function BCUT2D_MWHI>), ('BCUT2D_MWLOW', <function BCUT2D_MWLOW>), ('BCUT2D_CHGHI', <function BCUT2D_CHGHI>), ('BCUT2D_CHGLO', <function BCUT2D_CHGLO>), ('BCUT2D_LOGPHI', <function BCUT2D_LOGPHI>), ('BCUT2D_LOGPLOW', <function BCUT2D_LOGPLOW>), ('BCUT2D_MRHI', <function BCUT2D_MRHI>), ('BCUT2D_MRLOW', <function BCUT2D_MRLOW>), ('AvgIpc', <function AvgIpc>), ('BalabanJ', <function BalabanJ>), ('BertzCT', <function BertzCT>), ('Chi0', <function Chi0>), ('Chi0n', <function <lambda>>), ('Chi0v', <function <lambda>>), ('Chi1', <function Chi1>), ('Chi1n', <function <lambda>>), ('Chi1v', <function <lambda>>), ('Chi2n', <function <lambda>>), ('Chi2v', <function <lambda>>), ('Chi3n', <function <lambda>>), ('Chi3v', <function <lambda>>), ('Chi4n', <function <lambda>>), ('Chi4v', <function <lambda>>), ('HallKierAlpha', <function <lambda>>), ('Ipc', <function Ipc>), ('Kappa1', <function <lambda>>), ('Kappa2', <function <lambda>>), ('Kappa3', <function <lambda>>), ('LabuteASA', <function <lambda>>), ('PEOE_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA11', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA12', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA13', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA14', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA11', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA12', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('TPSA', <function <lambda>>), ('EState_VSA1', <function EState_VSA1>), ('EState_VSA10', <function EState_VSA10>), ('EState_VSA11', <function EState_VSA11>), ('EState_VSA2', <function EState_VSA2>), ('EState_VSA3', <function EState_VSA3>), ('EState_VSA4', <function EState_VSA4>), ('EState_VSA5', <function EState_VSA5>), ('EState_VSA6', <function EState_VSA6>), ('EState_VSA7', <function EState_VSA7>), ('EState_VSA8', <function EState_VSA8>), ('EState_VSA9', <function EState_VSA9>), ('VSA_EState1', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState10', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState2', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState3', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState4', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState5', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState6', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState7', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState8', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState9', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('FractionCSP3', <function <lambda>>), ('HeavyAtomCount', <function HeavyAtomCount>), ('NHOHCount', <function <lambda>>), ('NOCount', <function <lambda>>), ('NumAliphaticCarbocycles', <function <lambda>>), ('NumAliphaticHeterocycles', <function <lambda>>), ('NumAliphaticRings', <function <lambda>>), ('NumAromaticCarbocycles', <function <lambda>>), ('NumAromaticHeterocycles', <function <lambda>>), ('NumAromaticRings', <function <lambda>>), ('NumHAcceptors', <function <lambda>>), ('NumHDonors', <function <lambda>>), ('NumHeteroatoms', <function <lambda>>), ('NumRotatableBonds', <function <lambda>>), ('NumSaturatedCarbocycles', <function <lambda>>), ('NumSaturatedHeterocycles', <function <lambda>>), ('NumSaturatedRings', <function <lambda>>), ('RingCount', <function <lambda>>), ('MolLogP', <function <lambda>>), ('MolMR', <function <lambda>>), ('fr_Al_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH_noTert', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ArN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_N', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_NH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_O_noCOO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_S', <function _LoadPatterns.<locals>.<lambda>>), ('fr_HOCCN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Imine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH0', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_N_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Nhpyrrole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_SH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aldehyde', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_carbamate', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_halide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_allylic_oxid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amidine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aniline', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aryl_methyl', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_barbitur', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzodiazepine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_bicyclic', <function _LoadPatterns.<locals>.<lambda>>), ('fr_diazo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_dihydropyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_epoxide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ether', <function _LoadPatterns.<locals>.<lambda>>), ('fr_furan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_guanido', <function _LoadPatterns.<locals>.<lambda>>), ('fr_halogen', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imidazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isothiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone_Topliss', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactam', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_methoxy', <function _LoadPatterns.<locals>.<lambda>>), ('fr_morpholine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitrile', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom_nonortho', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitroso', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxime', <function _LoadPatterns.<locals>.<lambda>>), ('fr_para_hydroxylation', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol_noOrthoHbond', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_acid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperdine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_priamide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_prisulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_pyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_quatN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_term_acetylene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_tetrazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiophene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_unbrch_alkane', <function _LoadPatterns.<locals>.<lambda>>), ('fr_urea', <function _LoadPatterns.<locals>.<lambda>>)]
autocorr: rdkit.Chem.ChemUtils.DescriptorUtilities.VectorDescriptorWrapper
descList: list # value = [('MaxAbsEStateIndex', <function MaxAbsEStateIndex>), ('MaxEStateIndex', <function MaxEStateIndex>), ('MinAbsEStateIndex', <function MinAbsEStateIndex>), ('MinEStateIndex', <function MinEStateIndex>), ('qed', <function qed>), ('MolWt', <function <lambda>>), ('HeavyAtomMolWt', <function HeavyAtomMolWt>), ('ExactMolWt', <function <lambda>>), ('NumValenceElectrons', <function NumValenceElectrons>), ('NumRadicalElectrons', <function NumRadicalElectrons>), ('MaxPartialCharge', <function MaxPartialCharge>), ('MinPartialCharge', <function MinPartialCharge>), ('MaxAbsPartialCharge', <function MaxAbsPartialCharge>), ('MinAbsPartialCharge', <function MinAbsPartialCharge>), ('FpDensityMorgan1', <function FpDensityMorgan1>), ('FpDensityMorgan2', <function FpDensityMorgan2>), ('FpDensityMorgan3', <function FpDensityMorgan3>), ('BCUT2D_MWHI', <function BCUT2D_MWHI>), ('BCUT2D_MWLOW', <function BCUT2D_MWLOW>), ('BCUT2D_CHGHI', <function BCUT2D_CHGHI>), ('BCUT2D_CHGLO', <function BCUT2D_CHGLO>), ('BCUT2D_LOGPHI', <function BCUT2D_LOGPHI>), ('BCUT2D_LOGPLOW', <function BCUT2D_LOGPLOW>), ('BCUT2D_MRHI', <function BCUT2D_MRHI>), ('BCUT2D_MRLOW', <function BCUT2D_MRLOW>), ('AvgIpc', <function AvgIpc>), ('BalabanJ', <function BalabanJ>), ('BertzCT', <function BertzCT>), ('Chi0', <function Chi0>), ('Chi0n', <function <lambda>>), ('Chi0v', <function <lambda>>), ('Chi1', <function Chi1>), ('Chi1n', <function <lambda>>), ('Chi1v', <function <lambda>>), ('Chi2n', <function <lambda>>), ('Chi2v', <function <lambda>>), ('Chi3n', <function <lambda>>), ('Chi3v', <function <lambda>>), ('Chi4n', <function <lambda>>), ('Chi4v', <function <lambda>>), ('HallKierAlpha', <function <lambda>>), ('Ipc', <function Ipc>), ('Kappa1', <function <lambda>>), ('Kappa2', <function <lambda>>), ('Kappa3', <function <lambda>>), ('LabuteASA', <function <lambda>>), ('PEOE_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA11', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA12', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA13', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA14', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('PEOE_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('SMR_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA1', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA10', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA11', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA12', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA2', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA3', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA4', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA5', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA6', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA7', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA8', <function _InstallDescriptors.<locals>.<lambda>>), ('SlogP_VSA9', <function _InstallDescriptors.<locals>.<lambda>>), ('TPSA', <function <lambda>>), ('EState_VSA1', <function EState_VSA1>), ('EState_VSA10', <function EState_VSA10>), ('EState_VSA11', <function EState_VSA11>), ('EState_VSA2', <function EState_VSA2>), ('EState_VSA3', <function EState_VSA3>), ('EState_VSA4', <function EState_VSA4>), ('EState_VSA5', <function EState_VSA5>), ('EState_VSA6', <function EState_VSA6>), ('EState_VSA7', <function EState_VSA7>), ('EState_VSA8', <function EState_VSA8>), ('EState_VSA9', <function EState_VSA9>), ('VSA_EState1', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState10', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState2', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState3', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState4', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState5', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState6', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState7', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState8', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('VSA_EState9', <function _descriptor_VSA_EState.<locals>.VSA_EState_bin>), ('FractionCSP3', <function <lambda>>), ('HeavyAtomCount', <function HeavyAtomCount>), ('NHOHCount', <function <lambda>>), ('NOCount', <function <lambda>>), ('NumAliphaticCarbocycles', <function <lambda>>), ('NumAliphaticHeterocycles', <function <lambda>>), ('NumAliphaticRings', <function <lambda>>), ('NumAromaticCarbocycles', <function <lambda>>), ('NumAromaticHeterocycles', <function <lambda>>), ('NumAromaticRings', <function <lambda>>), ('NumHAcceptors', <function <lambda>>), ('NumHDonors', <function <lambda>>), ('NumHeteroatoms', <function <lambda>>), ('NumRotatableBonds', <function <lambda>>), ('NumSaturatedCarbocycles', <function <lambda>>), ('NumSaturatedHeterocycles', <function <lambda>>), ('NumSaturatedRings', <function <lambda>>), ('RingCount', <function <lambda>>), ('MolLogP', <function <lambda>>), ('MolMR', <function <lambda>>), ('fr_Al_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Al_OH_noTert', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ArN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_N', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_NH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ar_OH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_COO2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_O_noCOO', <function _LoadPatterns.<locals>.<lambda>>), ('fr_C_S', <function _LoadPatterns.<locals>.<lambda>>), ('fr_HOCCN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Imine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH0', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_NH2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_N_O', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation1', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Ndealkylation2', <function _LoadPatterns.<locals>.<lambda>>), ('fr_Nhpyrrole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_SH', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aldehyde', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_carbamate', <function _LoadPatterns.<locals>.<lambda>>), ('fr_alkyl_halide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_allylic_oxid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_amidine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aniline', <function _LoadPatterns.<locals>.<lambda>>), ('fr_aryl_methyl', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_azo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_barbitur', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_benzodiazepine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_bicyclic', <function _LoadPatterns.<locals>.<lambda>>), ('fr_diazo', <function _LoadPatterns.<locals>.<lambda>>), ('fr_dihydropyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_epoxide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ether', <function _LoadPatterns.<locals>.<lambda>>), ('fr_furan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_guanido', <function _LoadPatterns.<locals>.<lambda>>), ('fr_halogen', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_hdrzone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imidazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_imide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_isothiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_ketone_Topliss', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactam', <function _LoadPatterns.<locals>.<lambda>>), ('fr_lactone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_methoxy', <function _LoadPatterns.<locals>.<lambda>>), ('fr_morpholine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitrile', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitro_arom_nonortho', <function _LoadPatterns.<locals>.<lambda>>), ('fr_nitroso', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_oxime', <function _LoadPatterns.<locals>.<lambda>>), ('fr_para_hydroxylation', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phenol_noOrthoHbond', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_acid', <function _LoadPatterns.<locals>.<lambda>>), ('fr_phos_ester', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperdine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_piperzine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_priamide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_prisulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_pyridine', <function _LoadPatterns.<locals>.<lambda>>), ('fr_quatN', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfide', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfonamd', <function _LoadPatterns.<locals>.<lambda>>), ('fr_sulfone', <function _LoadPatterns.<locals>.<lambda>>), ('fr_term_acetylene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_tetrazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiazole', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiocyan', <function _LoadPatterns.<locals>.<lambda>>), ('fr_thiophene', <function _LoadPatterns.<locals>.<lambda>>), ('fr_unbrch_alkane', <function _LoadPatterns.<locals>.<lambda>>), ('fr_urea', <function _LoadPatterns.<locals>.<lambda>>)]
names = ['AUTOCORR2D_1', 'AUTOCORR2D_2', 'AUTOCORR2D_3', 'AUTOCORR2D_4', 'AUTOCORR2D_5', 'AUTOCORR2D_6', 'AUTOCORR2D_7', 'AUTOCORR2D_8', 'AUTOCORR2D_9', 'AUTOCORR2D_10', 'AUTOCORR2D_11', 'AUTOCORR2D_12', 'AUTOCORR2D_13', 'AUTOCORR2D_14', 'AUTOCORR2D_15', 'AUTOCORR2D_16', 'AUTOCORR2D_17', 'AUTOCORR2D_18', 'AUTOCORR2D_19', 'AUTOCORR2D_20', 'AUTOCORR2D_21', 'AUTOCORR2D_22', 'AUTOCORR2D_23', 'AUTOCORR2D_24', 'AUTOCORR2D_25', 'AUTOCORR2D_26', 'AUTOCORR2D_27', 'AUTOCORR2D_28', 'AUTOCORR2D_29', 'AUTOCORR2D_30', 'AUTOCORR2D_31', 'AUTOCORR2D_32', 'AUTOCORR2D_33', 'AUTOCORR2D_34', 'AUTOCORR2D_35', 'AUTOCORR2D_36', 'AUTOCORR2D_37', 'AUTOCORR2D_38', 'AUTOCORR2D_39', 'AUTOCORR2D_40', 'AUTOCORR2D_41', 'AUTOCORR2D_42', 'AUTOCORR2D_43', 'AUTOCORR2D_44', 'AUTOCORR2D_45', 'AUTOCORR2D_46', 'AUTOCORR2D_47', 'AUTOCORR2D_48', 'AUTOCORR2D_49', 'AUTOCORR2D_50', 'AUTOCORR2D_51', 'AUTOCORR2D_52', 'AUTOCORR2D_53', 'AUTOCORR2D_54', 'AUTOCORR2D_55', 'AUTOCORR2D_56', 'AUTOCORR2D_57', 'AUTOCORR2D_58', 'AUTOCORR2D_59', 'AUTOCORR2D_60', 'AUTOCORR2D_61', 'AUTOCORR2D_62', 'AUTOCORR2D_63', 'AUTOCORR2D_64', 'AUTOCORR2D_65', 'AUTOCORR2D_66', 'AUTOCORR2D_67', 'AUTOCORR2D_68', 'AUTOCORR2D_69', 'AUTOCORR2D_70', 'AUTOCORR2D_71', 'AUTOCORR2D_72', 'AUTOCORR2D_73', 'AUTOCORR2D_74', 'AUTOCORR2D_75', 'AUTOCORR2D_76', 'AUTOCORR2D_77', 'AUTOCORR2D_78', 'AUTOCORR2D_79', 'AUTOCORR2D_80', 'AUTOCORR2D_81', 'AUTOCORR2D_82', 'AUTOCORR2D_83', 'AUTOCORR2D_84', 'AUTOCORR2D_85', 'AUTOCORR2D_86', 'AUTOCORR2D_87', 'AUTOCORR2D_88', 'AUTOCORR2D_89', 'AUTOCORR2D_90', 'AUTOCORR2D_91', 'AUTOCORR2D_92', 'AUTOCORR2D_93', 'AUTOCORR2D_94', 'AUTOCORR2D_95', 'AUTOCORR2D_96', 'AUTOCORR2D_97', 'AUTOCORR2D_98', 'AUTOCORR2D_99', 'AUTOCORR2D_100', 'AUTOCORR2D_101', 'AUTOCORR2D_102', 'AUTOCORR2D_103', 'AUTOCORR2D_104', 'AUTOCORR2D_105', 'AUTOCORR2D_106', 'AUTOCORR2D_107', 'AUTOCORR2D_108', 'AUTOCORR2D_109', 'AUTOCORR2D_110', 'AUTOCORR2D_111', 'AUTOCORR2D_112', 'AUTOCORR2D_113', 'AUTOCORR2D_114', 'AUTOCORR2D_115', 'AUTOCORR2D_116', 'AUTOCORR2D_117', 'AUTOCORR2D_118', 'AUTOCORR2D_119', 'AUTOCORR2D_120', 'AUTOCORR2D_121', 'AUTOCORR2D_122', 'AUTOCORR2D_123', 'AUTOCORR2D_124', 'AUTOCORR2D_125', 'AUTOCORR2D_126', 'AUTOCORR2D_127', 'AUTOCORR2D_128', 'AUTOCORR2D_129', 'AUTOCORR2D_130', 'AUTOCORR2D_131', 'AUTOCORR2D_132', 'AUTOCORR2D_133', 'AUTOCORR2D_134', 'AUTOCORR2D_135', 'AUTOCORR2D_136', 'AUTOCORR2D_137', 'AUTOCORR2D_138', 'AUTOCORR2D_139', 'AUTOCORR2D_140', 'AUTOCORR2D_141', 'AUTOCORR2D_142', 'AUTOCORR2D_143', 'AUTOCORR2D_144', 'AUTOCORR2D_145', 'AUTOCORR2D_146', 'AUTOCORR2D_147', 'AUTOCORR2D_148', 'AUTOCORR2D_149', 'AUTOCORR2D_150', 'AUTOCORR2D_151', 'AUTOCORR2D_152', 'AUTOCORR2D_153', 'AUTOCORR2D_154', 'AUTOCORR2D_155', 'AUTOCORR2D_156', 'AUTOCORR2D_157', 'AUTOCORR2D_158', 'AUTOCORR2D_159', 'AUTOCORR2D_160', 'AUTOCORR2D_161', 'AUTOCORR2D_162', 'AUTOCORR2D_163', 'AUTOCORR2D_164', 'AUTOCORR2D_165', 'AUTOCORR2D_166', 'AUTOCORR2D_167', 'AUTOCORR2D_168', 'AUTOCORR2D_169', 'AUTOCORR2D_170', 'AUTOCORR2D_171', 'AUTOCORR2D_172', 'AUTOCORR2D_173', 'AUTOCORR2D_174', 'AUTOCORR2D_175', 'AUTOCORR2D_176', 'AUTOCORR2D_177', 'AUTOCORR2D_178', 'AUTOCORR2D_179', 'AUTOCORR2D_180', 'AUTOCORR2D_181', 'AUTOCORR2D_182', 'AUTOCORR2D_183', 'AUTOCORR2D_184', 'AUTOCORR2D_185', 'AUTOCORR2D_186', 'AUTOCORR2D_187', 'AUTOCORR2D_188', 'AUTOCORR2D_189', 'AUTOCORR2D_190', 'AUTOCORR2D_191', 'AUTOCORR2D_192']
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
ExactMolWt = rdkit.Chem.Descriptors.<lambda>
FractionCSP3 = rdkit.Chem.Lipinski.<lambda>
HallKierAlpha = rdkit.Chem.GraphDescriptors.<lambda>
Kappa1 = rdkit.Chem.GraphDescriptors.<lambda>
Kappa2 = rdkit.Chem.GraphDescriptors.<lambda>
Kappa3 = rdkit.Chem.GraphDescriptors.<lambda>
LabuteASA = rdkit.Chem.MolSurf.<lambda>
MolLogP = rdkit.Chem.Crippen.<lambda>
MolMR = rdkit.Chem.Crippen.<lambda>
MolWt = rdkit.Chem.Descriptors.<lambda>
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
PEOE_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA11 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA12 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA13 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA14 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
PEOE_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
RingCount = rdkit.Chem.Lipinski.<lambda>
SMR_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SMR_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA1 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA10 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA11 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA12 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA2 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA3 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA4 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA5 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA6 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA7 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA8 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
SlogP_VSA9 = rdkit.Chem.MolSurf._InstallDescriptors.<locals>.<lambda>
TPSA = rdkit.Chem.MolSurf.<lambda>
VSA_EState1 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState10 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState2 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState3 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState4 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState5 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState6 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState7 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState8 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
VSA_EState9 = rdkit.Chem.EState.EState_VSA._descriptor_VSA_EState.<locals>.VSA_EState_bin
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
