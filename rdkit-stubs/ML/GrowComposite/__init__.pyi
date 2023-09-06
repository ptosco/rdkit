""" command line utility for growing composite models

**Usage**

  _GrowComposite [optional args] filename_

**Command Line Arguments**

  - -n *count*: number of new models to build

  - -C *pickle file name*:  name of file containing composite upon which to build.

  - --inNote *note*: note to be used in loading composite models from the database
      for growing

  - --balTable *table name*:  table from which to take the original data set
     (for balancing)

  - --balWeight *weight*: (between 0 and 1) weighting factor for the new data
     (for balancing). OR, *weight* can be a list of weights

  - --balCnt *count*: number of individual models in the balanced composite
     (for balancing)

  - --balH: use only the holdout set from the original data set in the balancing
     (for balancing)

  - --balT: use only the training set from the original data set in the balancing
     (for balancing)

  - -S: shuffle the original data set
     (for balancing)

  - -r: randomize the activities of the original data set
     (for balancing)

  - -N *note*: note to be attached to the grown composite when it's saved in the
     database

  - --outNote *note*: equivalent to -N

  - -o *filename*: name of an output file to hold the pickled composite after
     it has been grown.
     If multiple balance weights are used, the weights will be added to
     the filenames.

  - -L *limit*: provide an (integer) limit on individual model complexity

  - -d *database name*: instead of reading the data from a QDAT file,
     pull it from a database.  In this case, the _filename_ argument
     provides the name of the database table containing the data set.

  - -p *tablename*: store persistence data in the database
     in table *tablename*

  - -l: locks the random number generator to give consistent sets
     of training and hold-out data.  This is primarily intended
     for testing purposes.

  - -g: be less greedy when training the models.

  - -G *number*: force trees to be rooted at descriptor *number*.

  - -D: show a detailed breakdown of the composite model performance
     across the training and, when appropriate, hold-out sets.

  - -t *threshold value*: use high-confidence predictions for the final
     analysis of the hold-out data.

  - -q *list string*:  Add QuantTrees to the composite and use the list
     specified in *list string* as the number of target quantization
     bounds for each descriptor.  Don't forget to include 0's at the
     beginning and end of *list string* for the name and value fields.
     For example, if there are 4 descriptors and you want 2 quant bounds
     apiece, you would use _-q "[0,2,2,2,2,0]"_.
     Two special cases:
       1) If you would like to ignore a descriptor in the model building,
          use '-1' for its number of quant bounds.
       2) If you have integer valued data that should not be quantized
          further, enter 0 for that descriptor.

  - -V: print the version number and exit

"""
from __future__ import annotations
import rdkit.ML.GrowComposite
import typing
from rdkit.Dbase.DbConnection import DbConnect
import numpy
import pickle
import rdkit.ML.BuildComposite
import rdkit.ML.Composite.AdjustComposite
import rdkit.ML.CompositeRun
import rdkit.ML.Data.DataUtils
import rdkit.ML.Data.SplitData
import rdkit.ML.ScreenComposite
import sys
import time
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AdjustComposite",
    "BalanceComposite",
    "BuildComposite",
    "CompositeRun",
    "DataUtils",
    "DbConnect",
    "GetComposites",
    "GrowIt",
    "ParseArgs",
    "ScreenComposite",
    "SetDefaults",
    "ShowVersion",
    "SplitData",
    "Usage",
    "message",
    "numpy",
    "pickle",
    "sys",
    "time"
]


__VERSION_STRING = '0.5.0'
_runDetails: rdkit.ML.CompositeRun.CompositeRun
_verbose = 1
