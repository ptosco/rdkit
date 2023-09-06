from __future__ import annotations
import rdkit.SimDivFilters.SimilarityPickers
import typing
from rdkit.DataStructs.TopNContainer import TopNContainer
import bisect
import rdkit.DataStructs

__all__ = [
    "DataStructs",
    "GenericPicker",
    "SpreadPicker",
    "TopNContainer",
    "TopNOverallPicker",
    "bisect"
]


class GenericPicker():
    _picks = None
    pass
class SpreadPicker(GenericPicker):
    """
    A class for picking the best matches across a library

    Connect to a database:

    >>> from rdkit import Chem
    >>> from rdkit import RDConfig
    >>> import os.path
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(dbName,'simple_mols1')
    >>> [x.upper() for x in conn.GetColumnNames()]
    ['SMILES', 'ID']
    >>> mols = []
    >>> for smi,id in conn.GetData():
    ...   mol = Chem.MolFromSmiles(str(smi))
    ...   mol.SetProp('_Name',str(id))
    ...   mols.append(mol)
    >>> len(mols)
    12

    Calculate fingerprints:

    >>> probefps = []
    >>> for mol in mols:
    ...   fp = Chem.RDKFingerprint(mol)
    ...   fp._id = mol.GetProp('_Name')
    ...   probefps.append(fp)

    Start by finding the top matches for a single probe.  This ether should pull
    other ethers from the db:

    >>> mol = Chem.MolFromSmiles('COC')
    >>> probeFp = Chem.RDKFingerprint(mol)
    >>> picker = SpreadPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
    >>> len(picker)
    2
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0

    The results come back in order:

    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    Now find the top matches for 2 probes.  We'll get one ether and one acid:

    >>> fps = []
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
    >>> picker = SpreadPicker(numToPick=3,probeFps=fps,dataSet=probefps)
    >>> len(picker)
    3
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'acid-1'
    >>> score
    1.0
    >>> fp,score = picker[2]
    >>> id = fp._id
    >>> str(id)
    'ether-2'
    """
    pass
class TopNOverallPicker(GenericPicker):
    """
    A class for picking the top N overall best matches across a library

    Connect to a database and build molecules:

    >>> from rdkit import Chem
    >>> from rdkit import RDConfig
    >>> import os.path
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(dbName,'simple_mols1')
    >>> [x.upper() for x in conn.GetColumnNames()]
    ['SMILES', 'ID']
    >>> mols = []
    >>> for smi,id in conn.GetData():
    ...   mol = Chem.MolFromSmiles(str(smi))
    ...   mol.SetProp('_Name',str(id))
    ...   mols.append(mol)
    >>> len(mols)
    12

    Calculate fingerprints:

    >>> probefps = []
    >>> for mol in mols:
    ...   fp = Chem.RDKFingerprint(mol)
    ...   fp._id = mol.GetProp('_Name')
    ...   probefps.append(fp)

    Start by finding the top matches for a single probe.  This ether should pull
    other ethers from the db:

    >>> mol = Chem.MolFromSmiles('COC')
    >>> probeFp = Chem.RDKFingerprint(mol)
    >>> picker = TopNOverallPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
    >>> len(picker)
    2
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0

    The results come back in order:

    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-2'

    Now find the top matches for 2 probes.  We'll get one ether and one acid:

    >>> fps = []
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
    >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
    >>> picker = TopNOverallPicker(numToPick=3,probeFps=fps,dataSet=probefps)
    >>> len(picker)
    3
    >>> fp,score = picker[0]
    >>> id = fp._id
    >>> str(id)
    'acid-1'
    >>> fp,score = picker[1]
    >>> id = fp._id
    >>> str(id)
    'ether-1'
    >>> score
    1.0
    >>> fp,score = picker[2]
    >>> id = fp._id
    >>> str(id)
    'acid-2'
    """
    pass
