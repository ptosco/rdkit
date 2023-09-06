""" descriptor calculator for compounds defined by a composition alone
  (only the composition is required)

"""
from __future__ import annotations
import rdkit.ML.Descriptors.CompoundDescriptors
import typing
import rdkit.ML.Descriptors.Descriptors
import rdkit.ML.Descriptors.Parser
import rdkit.RDConfig
import rdkit.utils.chemutils

__all__ = [
    "CompoundDescriptorCalculator",
    "Descriptors",
    "GetAllDescriptorNames",
    "Parser",
    "RDConfig",
    "chemutils",
    "countOptions"
]


class CompoundDescriptorCalculator(rdkit.ML.Descriptors.Descriptors.DescriptorCalculator):
    """
    used for calculating descriptors

      This is the central point for descriptor calculation

      **Notes**

      - There are two kinds of descriptors this cares about:

         1) *Simple Descriptors* can be calculated solely using atomic descriptor
            values and the composition of the compound.  The full list of possible
            simple descriptors is determined by the types of *Calculator Methods*
            (see below) and the contents of an atomic database.

            Simple Descriptors can be marked as *nonZeroDescriptors*.  These are used
            to winnow out atom types where particular atomic descriptors are zero
            (usually indicating that the value is unknown)

            Simple Descriptors are maintained locally in the _simpleList_

         2) *Compound Descriptors* may rely upon more complicated computation schemes
            and descriptors for the compound as a whole (e.g. structural variables, etc.).
            The full list of compound descriptors is limitless.  They are calculated using
            the _ML.Descriptors.Parser_ module.

            Compound Descriptors are maintained locally in the _compoundList_

      - This class has a some special methods which are labelled as *Calculator Method*
        These are used internally to take atomic descriptors and reduce them to a single
        simple descriptor value for a composition.  They are primarily intended for internal use.

      - a *composition vector* is a list of 2-tuples: '[(atom1name,atom1Num),...]'
        where atom1Num is the contribution of the atom to the stoichiometry of the
        compound. No assumption is made about the stoichiometries (i.e. they don't
        have to be either integral or all sum to one).

     
    """
    pass
countOptions = [('NVAL', 'total number of valence electrons'), ('NVAL_NO_FULL_F', 'number of valence electrons neglecting filled f shells'), ('NVAL_NO_FULL_D', 'number of valence electrons neglecting filled d shells'), ('NVAL_NO_FULL', 'number of valence electrons neglecting filled f and d shells')]
