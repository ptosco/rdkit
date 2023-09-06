""" The "parser" for compound descriptors.

I almost hesitate to document this, because it's not the prettiest
thing the world has ever seen... but it does work (for at least some
definitions of the word).

Rather than getting into the whole mess of writing a parser for the
compound descriptor expressions, I'm just using string substitutions
and python's wonderful ability to *eval* code.

It would probably be a good idea at some point to replace this with a
real parser, if only for the flexibility and intelligent error
messages that would become possible.

The general idea is that we're going to deal with expressions where
atomic descriptors have some kind of method applied to them which
reduces them to a single number for the entire composition.  Compound
descriptors (those applicable to the compound as a whole) are not
operated on by anything in particular (except for standard math stuff).

Here's the general flow of things:

  1) Composition descriptor references ($a, $b, etc.) are replaced with the
     corresponding descriptor names using string substitution.
     (*_SubForCompoundDescriptors*)

  2) Atomic descriptor references ($1, $2, etc) are replaced with lookups
     into the atomic dict with "DEADBEEF" in place of the atom name.
     (*_SubForAtomicVars*)

  3) Calls to Calculator Functions are augmented with a reference to
     the composition and atomic dictionary
     (*_SubMethodArgs*)

**NOTE:**

  anytime we don't know the answer for a descriptor, rather than
  throwing a (completely incomprehensible) exception, we just return
  -666.  So bad descriptor values should stand out like sore thumbs.

"""
from __future__ import annotations
import rdkit.ML.Descriptors.Parser
import typing
import rdkit.RDConfig

__all__ = [
    "AVG",
    "CalcMultipleCompoundsDescriptor",
    "CalcSingleCompoundDescriptor",
    "DEV",
    "HAS",
    "MAX",
    "MEAN",
    "MIN",
    "RDConfig",
    "SUM",
    "acos",
    "acosh",
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "ceil",
    "comb",
    "copysign",
    "cos",
    "cosh",
    "degrees",
    "dist",
    "e",
    "erf",
    "erfc",
    "exp",
    "expm1",
    "fabs",
    "factorial",
    "floor",
    "fmod",
    "frexp",
    "fsum",
    "gamma",
    "gcd",
    "hypot",
    "inf",
    "isclose",
    "isfinite",
    "isinf",
    "isnan",
    "isqrt",
    "knownMethods",
    "ldexp",
    "lgamma",
    "log",
    "log10",
    "log1p",
    "log2",
    "modf",
    "nan",
    "perm",
    "pi",
    "pow",
    "prod",
    "radians",
    "remainder",
    "sin",
    "sinh",
    "sqrt",
    "tan",
    "tanh",
    "tau",
    "trunc"
]


def hypot(*coordinates) -> value:
    """
    Multidimensional Euclidean distance from the origin to a point.

    Roughly equivalent to:
        sqrt(sum(x**2 for x in coordinates))

    For a two dimensional point (x, y), gives the hypotenuse
    using the Pythagorean theorem:  sqrt(x*x + y*y).

    For example, the hypotenuse of a 3/4/5 right triangle is:

        >>> hypot(3.0, 4.0)
        5.0
    """
__DEBUG = False
e = 2.718281828459045
inf: float # value = inf
knownMethods = ['SUM', 'MIN', 'MAX', 'MEAN', 'AVG', 'DEV', 'HAS']
nan: float # value = nan
pi = 3.141592653589793
tau = 6.283185307179586
AVG = rdkit.ML.Descriptors.Parser.MEAN
