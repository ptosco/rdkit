from __future__ import annotations
import rdkit.Chem.BRICS
import typing
import copy
import rdkit.Chem
import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdchem
import rdkit.RDRandom
import re
import sys

__all__ = [
    "BRICSBuild",
    "BRICSDecompose",
    "BreakBRICSBonds",
    "Chem",
    "FindBRICSBonds",
    "Reactions",
    "bType",
    "bnd",
    "bondMatchers",
    "compats",
    "copy",
    "defn",
    "dummyPattern",
    "e1",
    "e2",
    "env",
    "environMatchers",
    "environs",
    "g1",
    "g2",
    "gp",
    "i",
    "i1",
    "i2",
    "j",
    "labels",
    "patt",
    "ps",
    "r1",
    "r2",
    "random",
    "re",
    "reactionDefs",
    "reactions",
    "reverseReactions",
    "rs",
    "rxn",
    "rxnSet",
    "sma",
    "smartsGps",
    "sys",
    "t",
    "tmp"
]


bType = '-'
bnd = '-'
bondMatchers: list # value = [[('1', '3', '-', <rdkit.Chem.rdchem.Mol object>), ('1', '5', '-', <rdkit.Chem.rdchem.Mol object>), ('1', '10', '-', <rdkit.Chem.rdchem.Mol object>)], [('3', '4', '-', <rdkit.Chem.rdchem.Mol object>), ('3', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('3', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('3', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('3', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('4', '5', '-', <rdkit.Chem.rdchem.Mol object>), ('4', '11', '-', <rdkit.Chem.rdchem.Mol object>)], [('5', '12', '-', <rdkit.Chem.rdchem.Mol object>), ('5', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('5', '16', '-', <rdkit.Chem.rdchem.Mol object>), ('5', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('5', '15', '-', <rdkit.Chem.rdchem.Mol object>)], [('6', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('6', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('6', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('6', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('7a', '7b', '=', <rdkit.Chem.rdchem.Mol object>)], [('8', '9', '-', <rdkit.Chem.rdchem.Mol object>), ('8', '10', '-', <rdkit.Chem.rdchem.Mol object>), ('8', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('8', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('8', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('8', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('9', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('9', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('9', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('9', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('10', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('10', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('10', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('10', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('11', '13', '-', <rdkit.Chem.rdchem.Mol object>), ('11', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('11', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('11', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('13', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('13', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('13', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('14', '14', '-', <rdkit.Chem.rdchem.Mol object>), ('14', '15', '-', <rdkit.Chem.rdchem.Mol object>), ('14', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('15', '16', '-', <rdkit.Chem.rdchem.Mol object>)], [('16', '16', '-', <rdkit.Chem.rdchem.Mol object>)]]
compats = [('16', '16', '-')]
defn = '[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:2]'
dummyPattern: rdkit.Chem.rdchem.Mol
e1 = '[c;$(c(:c):c)]'
e2 = '[c;$(c(:c):c)]'
env = 'L16b'
environMatchers: dict # value = {'L1': <rdkit.Chem.rdchem.Mol object>, 'L3': <rdkit.Chem.rdchem.Mol object>, 'L4': <rdkit.Chem.rdchem.Mol object>, 'L5': <rdkit.Chem.rdchem.Mol object>, 'L6': <rdkit.Chem.rdchem.Mol object>, 'L7a': <rdkit.Chem.rdchem.Mol object>, 'L7b': <rdkit.Chem.rdchem.Mol object>, '#L8': <rdkit.Chem.rdchem.Mol object>, 'L8': <rdkit.Chem.rdchem.Mol object>, 'L9': <rdkit.Chem.rdchem.Mol object>, 'L10': <rdkit.Chem.rdchem.Mol object>, 'L11': <rdkit.Chem.rdchem.Mol object>, 'L12': <rdkit.Chem.rdchem.Mol object>, 'L13': <rdkit.Chem.rdchem.Mol object>, 'L14': <rdkit.Chem.rdchem.Mol object>, 'L14b': <rdkit.Chem.rdchem.Mol object>, 'L15': <rdkit.Chem.rdchem.Mol object>, 'L16': <rdkit.Chem.rdchem.Mol object>, 'L16b': <rdkit.Chem.rdchem.Mol object>}
environs = {'L1': '[C;D3]([#0,#6,#7,#8])(=O)', 'L3': '[O;D2]-;!@[#0,#6,#1]', 'L4': '[C;!D1;!$(C=*)]-;!@[#6]', 'L5': '[N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]', 'L6': '[C;D3;!R](=O)-;!@[#0,#6,#7,#8]', 'L7a': '[C;D2,D3]-[#6]', 'L7b': '[C;D2,D3]-[#6]', '#L8': '[C;!R;!D1]-;!@[#6]', 'L8': '[C;!R;!D1;!$(C!-*)]', 'L9': '[n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]', 'L10': '[N;R;$(N(@C(=O))@[C,N,O,S])]', 'L11': '[S;D2](-;!@[#0,#6])', 'L12': '[S;D4]([#6,#0])(=O)(=O)', 'L13': '[C;$(C(-;@[C,N,O,S])-;@[N,O,S])]', 'L14': '[c;$(c(:[c,n,o,s]):[n,o,s])]', 'L14b': '[c;$(c(:[c,n,o,s]):[n,o,s])]', 'L15': '[C;$(C(-;@C)-;@C)]', 'L16': '[c;$(c(:c):c)]', 'L16b': '[c;$(c(:c):c)]'}
g1 = '16'
g2 = '16'
gp = ['[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:2]']
i = 13
i1 = '16'
i2 = '16'
j = 0
labels = ['16', '16']
patt: rdkit.Chem.rdchem.Mol
ps = '[16*]-[*:1].[16*]-[*:2]'
r1 = '[c;$(c(:c):c)]'
r2 = '[c;$(c(:c):c)]'
reactionDefs = ([('1', '3', '-'), ('1', '5', '-'), ('1', '10', '-')], [('3', '4', '-'), ('3', '13', '-'), ('3', '14', '-'), ('3', '15', '-'), ('3', '16', '-')], [('4', '5', '-'), ('4', '11', '-')], [('5', '12', '-'), ('5', '14', '-'), ('5', '16', '-'), ('5', '13', '-'), ('5', '15', '-')], [('6', '13', '-'), ('6', '14', '-'), ('6', '15', '-'), ('6', '16', '-')], [('7a', '7b', '=')], [('8', '9', '-'), ('8', '10', '-'), ('8', '13', '-'), ('8', '14', '-'), ('8', '15', '-'), ('8', '16', '-')], [('9', '13', '-'), ('9', '14', '-'), ('9', '15', '-'), ('9', '16', '-')], [('10', '13', '-'), ('10', '14', '-'), ('10', '15', '-'), ('10', '16', '-')], [('11', '13', '-'), ('11', '14', '-'), ('11', '15', '-'), ('11', '16', '-')], [('13', '14', '-'), ('13', '15', '-'), ('13', '16', '-')], [('14', '14', '-'), ('14', '15', '-'), ('14', '16', '-')], [('15', '16', '-')], [('16', '16', '-')])
reactions: tuple # value = ([<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>], [<rdkit.Chem.rdChemReactions.ChemicalReaction object>])
reverseReactions: list # value = [<rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>, <rdkit.Chem.rdChemReactions.ChemicalReaction object>]
rs = '[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]'
rxn: rdkit.Chem.rdChemReactions.ChemicalReaction
rxnSet = ['[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:2]']
sma = '[16*]-[*:1].[16*]-[*:2]>>[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]'
smartsGps = (['[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([O;D2]-;!@[#0,#6,#1]):2]>>[1*]-[*:1].[3*]-[*:2]', '[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):2]>>[1*]-[*:1].[5*]-[*:2]', '[$([C;D3]([#0,#6,#7,#8])(=O)):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>>[1*]-[*:1].[10*]-[*:2]'], ['[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;!D1;!$(C=*)]-;!@[#6]):2]>>[3*]-[*:1].[4*]-[*:2]', '[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[3*]-[*:1].[13*]-[*:2]', '[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[3*]-[*:1].[14*]-[*:2]', '[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[3*]-[*:1].[15*]-[*:2]', '[$([O;D2]-;!@[#0,#6,#1]):1]-;!@[$([c;$(c(:c):c)]):2]>>[3*]-[*:1].[16*]-[*:2]'], ['[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):2]>>[4*]-[*:1].[5*]-[*:2]', '[$([C;!D1;!$(C=*)]-;!@[#6]):1]-;!@[$([S;D2](-;!@[#0,#6])):2]>>[4*]-[*:1].[11*]-[*:2]'], ['[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([S;D4]([#6,#0])(=O)(=O)):2]>>[5*]-[*:1].[12*]-[*:2]', '[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[5*]-[*:1].[14*]-[*:2]', '[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[5*]-[*:1].[16*]-[*:2]', '[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[5*]-[*:1].[13*]-[*:2]', '[$([N;!D1;!$(N=*);!$(N-[!#6;!#16;!#0;!#1]);!$([N;R]@[C;R]=O)]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[5*]-[*:1].[15*]-[*:2]'], ['[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[6*]-[*:1].[13*]-[*:2]', '[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[6*]-[*:1].[14*]-[*:2]', '[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[6*]-[*:1].[15*]-[*:2]', '[$([C;D3;!R](=O)-;!@[#0,#6,#7,#8]):1]-;!@[$([c;$(c(:c):c)]):2]>>[6*]-[*:1].[16*]-[*:2]'], ['[$([C;D2,D3]-[#6]):1]=;!@[$([C;D2,D3]-[#6]):2]>>[7*]-[*:1].[7*]-[*:2]'], ['[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):2]>>[8*]-[*:1].[9*]-[*:2]', '[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([N;R;$(N(@C(=O))@[C,N,O,S])]):2]>>[8*]-[*:1].[10*]-[*:2]', '[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[8*]-[*:1].[13*]-[*:2]', '[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[8*]-[*:1].[14*]-[*:2]', '[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[8*]-[*:1].[15*]-[*:2]', '[$([C;!R;!D1;!$(C!-*)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[8*]-[*:1].[16*]-[*:2]'], ['[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[9*]-[*:1].[13*]-[*:2]', '[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[9*]-[*:1].[14*]-[*:2]', '[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[9*]-[*:1].[15*]-[*:2]', '[$([n;+0;$(n(:[c,n,o,s]):[c,n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[9*]-[*:1].[16*]-[*:2]'], ['[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[10*]-[*:1].[13*]-[*:2]', '[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[10*]-[*:1].[14*]-[*:2]', '[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[10*]-[*:1].[15*]-[*:2]', '[$([N;R;$(N(@C(=O))@[C,N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[10*]-[*:1].[16*]-[*:2]'], ['[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):2]>>[11*]-[*:1].[13*]-[*:2]', '[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[11*]-[*:1].[14*]-[*:2]', '[$([S;D2](-;!@[#0,#6])):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[11*]-[*:1].[15*]-[*:2]', '[$([S;D2](-;!@[#0,#6])):1]-;!@[$([c;$(c(:c):c)]):2]>>[11*]-[*:1].[16*]-[*:2]'], ['[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[13*]-[*:1].[14*]-[*:2]', '[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[13*]-[*:1].[15*]-[*:2]', '[$([C;$(C(-;@[C,N,O,S])-;@[N,O,S])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[13*]-[*:1].[16*]-[*:2]'], ['[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:[c,n,o,s]):[n,o,s])]):2]>>[14*]-[*:1].[14*]-[*:2]', '[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([C;$(C(-;@C)-;@C)]):2]>>[14*]-[*:1].[15*]-[*:2]', '[$([c;$(c(:[c,n,o,s]):[n,o,s])]):1]-;!@[$([c;$(c(:c):c)]):2]>>[14*]-[*:1].[16*]-[*:2]'], ['[$([C;$(C(-;@C)-;@C)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[15*]-[*:1].[16*]-[*:2]'], ['[$([c;$(c(:c):c)]):1]-;!@[$([c;$(c(:c):c)]):2]>>[16*]-[*:1].[16*]-[*:2]'])
t: rdkit.Chem.rdChemReactions.ChemicalReaction
tmp: list # value = [('16', '16', '-', <rdkit.Chem.rdchem.Mol object>)]
