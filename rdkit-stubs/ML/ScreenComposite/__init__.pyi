""" command line utility for screening composite models

**Usage**

  _ScreenComposite [optional args] modelfile(s) datafile_

Unless indicated otherwise (via command line arguments), _modelfile_ is
a file containing a pickled composite model and _filename_ is a QDAT file.

**Command Line Arguments**

  - -t *threshold value(s)*: use high-confidence predictions for the final
     analysis of the hold-out data.  The threshold value can be either a single
     float or a list/tuple of floats.  All thresholds should be between
     0.0 and 1.0

  - -D: do a detailed screen.

  - -d *database name*: instead of reading the data from a QDAT file,
     pull it from a database.  In this case, the _datafile_ argument
     provides the name of the database table containing the data set.

  - -N *note*: use all models from the database which have this note.
               The modelfile argument should contain the name of the table
               with the models.

  - -H: screen only the hold out set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -T: screen only the training set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -E: do a detailed Error analysis.  This shows each misclassified
     point and the number of times it was missed across all screened
     composites.  If the --enrich argument is also provided, only compounds
     that have true activity value equal to the enrichment value will be
     used.

  - --enrich *enrichVal*: target "active" value to be used in calculating
     enrichments.

  - -A: show All predictions.

  - -S: shuffle activity values before screening

  - -R: randomize activity values before screening

  - -F *filter frac*: filters the data before training to change the
     distribution of activity values in the training set.  *filter frac*
     is the fraction of the training set that should have the target value.
     **See note in BuildComposite help about data filtering**

  - -v *filter value*: filters the data before training to change the
     distribution of activity values in the training set. *filter value*
     is the target value to use in filtering.
     **See note in BuildComposite help about data filtering**

  - -V: be verbose when screening multiple models

  - -h: show this message and exit

  - --OOB: Do out an "out-of-bag" generalization error estimate.  This only
      makes sense when applied to the original data set.

  - --pickleCol *colId*: index of the column containing a pickled value
      (used primarily for cases where fingerprints are used as descriptors)

  *** Options for making Prediction (Hanneke) Plots ***

  - --predPlot=<fileName>: triggers the generation of a Hanneke plot and
      sets the name of the .txt file which will hold the output data.
      A Gnuplot control file, <fileName>.gnu, will also be generated.

  - --predActTable=<name> (optional):  name of the database table
      containing activity values.  If this is not provided, activities
      will be read from the same table containing the screening data

  - --predActCol=<name> (optional):  name of the activity column. If not
      provided, the name of the last column in the activity table will
      be used.

  - --predLogScale (optional):  If provided, the x axis of the
      prediction plot (the activity axis) will be plotted using a log
      scale

  - --predShow: launch a gnuplot instance and display the prediction
      plot (the plot will still be written to disk).

  *** The following options are likely obsolete ***

  - -P: read pickled data.  The datafile argument should contain
     a pickled data set. *relevant only to qdat files*

  - -q: data are not quantized (the composite should take care of
     quantization itself if it requires quantized data). *relevant only to
     qdat files*



"""
from __future__ import annotations
import rdkit.ML.ScreenComposite
import typing
from rdkit.Dbase.DbConnection import DbConnect
import PIL.Image
import PIL.ImageDraw
import numpy
import os
import pickle
import rdkit.DataStructs
import rdkit.Dbase.DbModule
import rdkit.ML.CompositeRun
import rdkit.ML.Data.DataUtils
import rdkit.ML.Data.SplitData
import sys
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CalcEnrichment",
    "CollectResults",
    "CompositeRun",
    "DataStructs",
    "DataUtils",
    "DbConnect",
    "DbModule",
    "DetailedScreen",
    "GetScreenImage",
    "Go",
    "Image",
    "ImageDraw",
    "MakePredPlot",
    "ParseArgs",
    "PrepareDataFromDetails",
    "ScreenFromDetails",
    "ScreenIt",
    "ScreenToHtml",
    "SetDefaults",
    "ShowVersion",
    "ShowVoteResults",
    "SplitData",
    "Usage",
    "error",
    "hasPil",
    "message",
    "numpy",
    "os",
    "pickle",
    "sys"
]


__VERSION_STRING = '3.3.0'
_details: rdkit.ML.CompositeRun.CompositeRun
hasPil = 1
