"""
Importing pandasTools enables several features that allow for using RDKit molecules as columns of a
Pandas dataframe.
If the dataframe is containing a molecule format in a column (e.g. smiles), like in this example:

>>> from rdkit.Chem import PandasTools
>>> import pandas as pd
>>> import os
>>> from rdkit import RDConfig
>>> antibiotics = pd.DataFrame(columns=['Name','Smiles'])
>>> antibiotics = pd.concat([antibiotics, pd.DataFrame.from_records([{'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',
...   'Name':'Penicilline G'}])], ignore_index=True) #Penicilline G
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O',
...   'Name':'Tetracycline'}])], ignore_index=True) #Tetracycline
>>> antibiotics = pd.concat([antibiotics,pd.DataFrame.from_records([{
...   'Smiles':'CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C',
...   'Name':'Ampicilline'}])], ignore_index=True) #Ampicilline
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles']
>>> print(antibiotics)
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
1   Tetracycline  CC1(C2CC3C(C(=O)C(=C(C3(C(=O)C2=C(C4=C1C=CC=C4...
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...

a new column can be created holding the respective RDKit molecule objects. The fingerprint can be
included to accelerate substructure searches on the dataframe.

>>> PandasTools.AddMoleculeColumnToFrame(antibiotics,'Smiles','Molecule',includeFingerprints=True)
>>> print([str(x) for x in  antibiotics.columns])
['Name', 'Smiles', 'Molecule']

A substructure filter can be applied on the dataframe using the RDKit molecule column,
because the ">=" operator has been modified to work as a substructure check.
Such the antibiotics containing the beta-lactam ring "C1C(=O)NC1" can be obtained by

>>> beta_lactam = Chem.MolFromSmiles('C1C(=O)NC1')
>>> beta_lactam_antibiotics = antibiotics[antibiotics['Molecule'] >= beta_lactam]
>>> print(beta_lactam_antibiotics[['Name','Smiles']])
            Name                                             Smiles
0  Penicilline G    CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C
2  Ampicilline  CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O...


It is also possible to load an SDF file can be load into a dataframe.

>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile,smilesName='SMILES',molColName='Molecule',
...            includeFingerprints=True)
>>> frame.info # doctest: +SKIP
<bound method DataFrame.info of <class 'pandas.core.frame.DataFrame'>
Int64Index: 200 entries, 0 to 199
Data columns:
AMW                       200  non-null values
CLOGP                     200  non-null values
CP                        200  non-null values
CR                        200  non-null values
DAYLIGHT.FPG              200  non-null values
DAYLIGHT_CLOGP            200  non-null values
FP                        200  non-null values
ID                        200  non-null values
ISM                       200  non-null values
LIPINSKI_VIOLATIONS       200  non-null values
NUM_HACCEPTORS            200  non-null values
NUM_HDONORS               200  non-null values
NUM_HETEROATOMS           200  non-null values
NUM_LIPINSKIHACCEPTORS    200  non-null values
NUM_LIPINSKIHDONORS       200  non-null values
NUM_RINGS                 200  non-null values
NUM_ROTATABLEBONDS        200  non-null values
P1                        30  non-null values
SMILES                    200  non-null values
Molecule                  200  non-null values
dtypes: object(20)>

The standard ForwardSDMolSupplier keywords are also available:

>>> sdfFile = os.path.join(RDConfig.RDDataDir,'NCI/first_200.props.sdf')
>>> frame = PandasTools.LoadSDF(sdfFile, smilesName='SMILES', molColName='Molecule',
...            includeFingerprints=True, removeHs=False, strictParsing=True)

Conversion to html is quite easy:

>>> htm = frame.to_html() # doctest:
...
>>> str(htm[:36])
'<table border="1" class="dataframe">'

In order to support rendering the molecules as images in the HTML export of the
dataframe, we use a custom formatter for columns containing RDKit molecules,
and also disable escaping of HTML where needed.
"""
from __future__ import annotations
import rdkit.Chem.PandasTools
import typing
from _io import BytesIO
from rdkit.Chem.rdmolfiles import SDWriter
import logging
import numpy
import pandas
import rdkit
import rdkit.Avalon.pyAvalonTools
import rdkit.Chem
import rdkit.Chem.AllChem
import rdkit.Chem.Draw
import rdkit.Chem.PandasPatcher
import rdkit.Chem.Scaffolds.MurckoScaffold
import rdkit.Chem.rdchem
import rdkit.DataStructs
import sys
import xml.dom.minidom
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AddMoleculeColumnToFrame",
    "AddMurckoToFrame",
    "AlignMol",
    "AlignToScaffold",
    "AllChem",
    "BytesIO",
    "ChangeMoleculeRendering",
    "Chem",
    "DataStructs",
    "Draw",
    "FrameToGridImage",
    "InstallPandasTools",
    "InteractiveRenderer",
    "LoadSDF",
    "MurckoScaffold",
    "PandasPatcher",
    "PrintAsImageString",
    "RGroupDecompositionToFrame",
    "RemoveSaltsFromFrame",
    "RenderImagesInAllDataFrames",
    "SDWriter",
    "SaveSMILESFromFrame",
    "SaveXlsxFromFrame",
    "UninstallPandasTools",
    "WriteSDF",
    "b64encode",
    "drawOptions",
    "highlightSubstructures",
    "log",
    "logging",
    "minidom",
    "molJustify",
    "molRepresentation",
    "molSize",
    "np",
    "pd",
    "pyAvalonTools",
    "rdchem",
    "rdkit",
    "sys"
]


InteractiveRenderer = None
_originalSettings: dict # value = {'Chem.Mol.__ge__': <slot wrapper '__ge__' of 'object' objects>}
_saltRemover = None
drawOptions = None
highlightSubstructures = True
log: logging.Logger # value = <Logger rdkit.Chem.PandasTools (WARNING)>
molJustify = 'center'
molRepresentation = 'png'
molSize = (200, 200)
