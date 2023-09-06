"""Command line tool to construct an enrichment plot from saved composite models

Usage:  EnrichPlot [optional args] -d dbname -t tablename <models>

Required Arguments:
  -d "dbName": the name of the database for screening

  -t "tablename": provide the name of the table with the data to be screened

  <models>: file name(s) of pickled composite model(s).
     If the -p argument is also provided (see below), this argument is ignored.

Optional Arguments:
  - -a "list": the list of result codes to be considered active.  This will be
        eval'ed, so be sure that it evaluates as a list or sequence of
        integers. For example, -a "[1,2]" will consider activity values 1 and 2
        to be active

  - --enrich "list": identical to the -a argument above.

  - --thresh: sets a threshold for the plot.  If the confidence falls below
          this value, picking will be terminated

  - -H: screen only the hold out set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

  - -T: screen only the training set (works only if a version of
        BuildComposite more recent than 1.2.2 was used).

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

  - -p "tableName": provides the name of a db table containing the
      models to be screened.  If you use this argument, you should also
      use the -N argument (below) to specify a note value.

  - -N "note": provides a note to be used to pull models from a db table.

  - --plotFile "filename": writes the data to an output text file (filename.dat)
    and creates a gnuplot input file (filename.gnu) to plot it

  - --showPlot: causes the gnuplot plot constructed using --plotFile to be
    displayed in gnuplot.

"""
from __future__ import annotations
import rdkit.ML.EnrichPlot
import typing
from rdkit.Dbase.DbConnection import DbConnect
import numpy
import pickle
import rdkit.DataStructs
import rdkit.ML.CompositeRun
import rdkit.ML.Data.DataUtils
import rdkit.ML.Data.SplitData
import rdkit.ML.Data.Stats
import rdkit.RDConfig
import sys
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AccumulateCounts",
    "CompositeRun",
    "DataStructs",
    "DataUtils",
    "DbConnect",
    "MakePlot",
    "RDConfig",
    "ScreenModel",
    "SplitData",
    "Stats",
    "Usage",
    "cmp",
    "error",
    "message",
    "numpy",
    "pickle",
    "sys"
]


__VERSION_STRING = '2.4.0'
