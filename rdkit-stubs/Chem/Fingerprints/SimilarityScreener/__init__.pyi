""" class definitions for similarity screening

See _SimilarityScreener_ for overview of required API

"""
from __future__ import annotations
import rdkit.Chem.Fingerprints.SimilarityScreener
import typing
import rdkit.DataStructs
import rdkit.DataStructs.TopNContainer

__all__ = [
    "DataStructs",
    "SimilarityScreener",
    "ThresholdScreener",
    "TopNContainer",
    "TopNScreener"
]


class SimilarityScreener():
    """
    base class

       important attributes:
          probe: the probe fingerprint against which we screen.

          metric: a function that takes two arguments and returns a similarity
                  measure between them

          dataSource: the source pool from which to draw, needs to support
                  a next() method

          fingerprinter: a function that takes a molecule and returns a
                 fingerprint of the appropriate format


        **Notes**
           subclasses must support either an iterator interface
           or __len__ and __getitem__
      
    """
    pass
class ThresholdScreener(SimilarityScreener):
    """
    Used to return all compounds that have a similarity
         to the probe beyond a threshold value

        **Notes**:

          - This is as lazy as possible, so the data source isn't
            queried until the client asks for a hit.

          - In addition to being lazy, this class is as thin as possible.
            (Who'd have thought it was possible!)
            Hits are *not* stored locally, so if a client resets
            the iteration and starts over, the same amount of work must
            be done to retrieve the hits.

          - The thinness and laziness forces us to support only forward
            iteration (not random access)

       
    """
    pass
class TopNScreener(SimilarityScreener):
    """
    A screener that only returns the top N hits found

         **Notes**

           - supports forward iteration and getitem

       
    """
    pass
