""" utility functionality for fingerprinting sets of molecules
 includes a command line app for working with fingerprints
 and databases


Sample Usage:

  python FingerprintMols.py  -d data.gdb         -t 'raw_dop_data' --smilesName="Structure" --idName="Mol_ID"          --outTable="daylight_sig"


"""
from __future__ import annotations
import rdkit.Chem.Fingerprints.FingerprintMols
import typing
import getopt
import pickle
import rdkit.Chem
import rdkit.Chem.MACCSkeys
import rdkit.DataStructs
import rdkit.ML.Cluster.Murtagh
import sys

__all__ = [
    "Chem",
    "DataStructs",
    "FingerprintMol",
    "FingerprinterDetails",
    "FingerprintsFromDetails",
    "FingerprintsFromMols",
    "FingerprintsFromPickles",
    "FingerprintsFromSmiles",
    "FoldFingerprintToTargetDensity",
    "GetRDKFingerprint",
    "MACCSkeys",
    "Murtagh",
    "ParseArgs",
    "Usage",
    "error",
    "getopt",
    "message",
    "pickle",
    "sys"
]


class FingerprinterDetails():
    """
    class for storing the details of a fingerprinting run,
        generates sensible defaults on construction

     
    """
    pass
_usageDoc = "\nUsage: FingerprintMols.py [args] <fName>\n\n  If <fName> is provided and no tableName is specified (see below),\n  data will be read from the text file <fName>.  Text files delimited\n  with either commas (extension .csv) or tabs (extension .txt) are\n  supported.\n\n  Command line arguments are:\n    - -d _dbName_: set the name of the database from which\n      to pull input molecule information.  If output is\n      going to a database, this will also be used for that\n      unless the --outDbName option is used.\n\n    - -t _tableName_: set the name of the database table\n      from which to pull input molecule information\n\n    - --smilesName=val: sets the name of the SMILES column\n      in the input database.  Default is *SMILES*.\n\n    - --useSD:  Assume that the input file is an SD file, not a SMILES\n       table.\n\n    - --idName=val: sets the name of the id column in the input\n      database.  Defaults to be the name of the first db column\n      (or *ID* for text files).\n\n    - -o _outFileName_:  name of the output file (output will\n      be a pickle file with one label,fingerprint entry for each\n      molecule).\n\n    - --outTable=val: name of the output db table used to store\n      fingerprints.  If this table already exists, it will be\n      replaced.\n\n    - --outDbName: name of output database, if it's being used.\n      Defaults to be the same as the input db.\n\n    - --fpColName=val: name to use for the column which stores\n      fingerprints (in pickled format) in the output db table.\n      Default is *AutoFragmentFP*\n\n    - --maxSize=val:  base size of the fingerprints to be generated\n      Default is *2048*\n\n    - --minSize=val: minimum size of the fingerprints to be generated\n      (limits the amount of folding that happens).  Default is *64*\n\n    - --density=val: target bit density in the fingerprint.  The\n      fingerprint will be folded until this density is\n      reached. Default is *0.3*\n\n    - --minPath=val:  minimum path length to be included in\n      fragment-based fingerprints. Default is *1*.\n\n    - --maxPath=val:  maximum path length to be included in\n      fragment-based fingerprints. Default is *7*.\n\n    - --nBitsPerHash: number of bits to be set in the output\n      fingerprint for each fragment. Default is *2*.\n\n    - --discrim: use of path-based discriminators to hash bits.\n      Default is *false*.\n\n    - -V: include valence information in the fingerprints\n      Default is *false*.\n\n    - -H: include Hs in the fingerprint\n      Default is *false*.\n\n    - --maxMols=val: sets the maximum number of molecules to be\n      fingerprinted.\n\n    - --useMACCS: use the public MACCS keys to do the fingerprinting\n      (instead of a daylight-type fingerprint)\n\n"
