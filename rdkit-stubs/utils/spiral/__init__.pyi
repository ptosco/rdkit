from __future__ import annotations
import rdkit.utils.spiral
import typing
from numpy import AxisError
from numpy import ComplexWarning
from numpy import DataSource
from numpy import ModuleDeprecationWarning
from numpy import RankWarning
from numpy import TooHardError
from numpy import VisibleDeprecationWarning
from numpy import bool_
from numpy import broadcast
from numpy import busdaycalendar
from numpy import bytes_
from numpy import character
from numpy import chararray
from numpy import complex128
from numpy import complex256
from numpy import complex64
from numpy import complexfloating
from numpy import datetime64
from numpy import dtype
from numpy import errstate
from numpy import finfo
from numpy import flatiter
from numpy import flexible
from numpy import float128
from numpy import float16
from numpy import float32
from numpy import float64
from numpy import floating
from numpy import format_parser
from numpy import generic
from numpy import iinfo
from numpy import inexact
from numpy import int16
from numpy import int32
from numpy import int64
from numpy import int8
from numpy import integer
from numpy import longlong
from numpy import matrix
from numpy import memmap
from numpy import ndarray
from numpy import ndenumerate
from numpy import ndindex
from numpy import nditer
from numpy import number
from numpy import object_
from numpy import poly1d
from numpy import recarray
from numpy import record
from numpy import signedinteger
from numpy import str_
from numpy import timedelta64
from numpy import ufunc
from numpy import uint16
from numpy import uint32
from numpy import uint64
from numpy import uint8
from numpy import ulonglong
from numpy import unsignedinteger
from numpy import vectorize
from numpy import void
import math
import numpy
import numpy.core.defchararray
import numpy.core.numerictypes
import numpy.core.records
import numpy.ctypeslib
import numpy.fft
import numpy.lib.index_tricks
import numpy.lib.scimath
import numpy.linalg
import numpy.ma
import numpy.random
import rdkit.sping.pid
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ALLOW_THREADS",
    "AxisError",
    "BUFSIZE",
    "CLIP",
    "ComplexWarning",
    "DataSource",
    "DrawSpiral",
    "ERR_CALL",
    "ERR_DEFAULT",
    "ERR_IGNORE",
    "ERR_LOG",
    "ERR_PRINT",
    "ERR_RAISE",
    "ERR_WARN",
    "FLOATING_POINT_SUPPORT",
    "FPE_DIVIDEBYZERO",
    "FPE_INVALID",
    "FPE_OVERFLOW",
    "FPE_UNDERFLOW",
    "False_",
    "Inf",
    "Infinity",
    "MAXDIMS",
    "MAY_SHARE_BOUNDS",
    "MAY_SHARE_EXACT",
    "ModuleDeprecationWarning",
    "NAN",
    "NINF",
    "NZERO",
    "NaN",
    "PINF",
    "PZERO",
    "RAISE",
    "RankWarning",
    "SHIFT_DIVIDEBYZERO",
    "SHIFT_INVALID",
    "SHIFT_OVERFLOW",
    "SHIFT_UNDERFLOW",
    "ScalarType",
    "TooHardError",
    "True_",
    "UFUNC_BUFSIZE_DEFAULT",
    "UFUNC_PYVALS_NAME",
    "VisibleDeprecationWarning",
    "WRAP",
    "absolute",
    "add",
    "add_docstring",
    "add_newdoc",
    "add_newdoc_ufunc",
    "all",
    "allclose",
    "alltrue",
    "amax",
    "amin",
    "angle",
    "any",
    "append",
    "apply_along_axis",
    "apply_over_axes",
    "arange",
    "arccos",
    "arccosh",
    "arcsin",
    "arcsinh",
    "arctan",
    "arctan2",
    "arctanh",
    "argmax",
    "argmin",
    "argpartition",
    "argsort",
    "argwhere",
    "around",
    "array",
    "array2string",
    "array_equal",
    "array_equiv",
    "array_repr",
    "array_split",
    "array_str",
    "asanyarray",
    "asarray",
    "asarray_chkfinite",
    "ascontiguousarray",
    "asfarray",
    "asfortranarray",
    "asmatrix",
    "atleast_1d",
    "atleast_2d",
    "atleast_3d",
    "average",
    "bartlett",
    "base_repr",
    "binary_repr",
    "bincount",
    "bitwise_and",
    "bitwise_not",
    "bitwise_or",
    "bitwise_xor",
    "blackman",
    "block",
    "bmat",
    "bool8",
    "bool_",
    "broadcast",
    "broadcast_arrays",
    "broadcast_shapes",
    "broadcast_to",
    "busday_count",
    "busday_offset",
    "busdaycalendar",
    "byte",
    "byte_bounds",
    "bytes0",
    "bytes_",
    "c_",
    "can_cast",
    "cast",
    "cbrt",
    "cdouble",
    "ceil",
    "cfloat",
    "char",
    "character",
    "chararray",
    "choose",
    "clip",
    "clongdouble",
    "clongfloat",
    "column_stack",
    "common_type",
    "compare_chararrays",
    "complex128",
    "complex256",
    "complex64",
    "complex_",
    "complexfloating",
    "compress",
    "concatenate",
    "conj",
    "conjugate",
    "convolve",
    "copy",
    "copysign",
    "copyto",
    "corrcoef",
    "correlate",
    "cos",
    "cosh",
    "count_nonzero",
    "cov",
    "cross",
    "csingle",
    "ctypeslib",
    "cumprod",
    "cumproduct",
    "cumsum",
    "datetime64",
    "datetime_as_string",
    "datetime_data",
    "deg2rad",
    "degrees",
    "delete",
    "deprecate",
    "deprecate_with_doc",
    "diag",
    "diag_indices",
    "diag_indices_from",
    "diagflat",
    "diagonal",
    "diff",
    "digitize",
    "disp",
    "divide",
    "divmod",
    "dot",
    "double",
    "dsplit",
    "dstack",
    "dtype",
    "e",
    "ediff1d",
    "einsum",
    "einsum_path",
    "emath",
    "empty",
    "empty_like",
    "equal",
    "errstate",
    "euler_gamma",
    "exp",
    "exp2",
    "expand_dims",
    "expm1",
    "extract",
    "eye",
    "fabs",
    "fastCopyAndTranspose",
    "fft",
    "fill_diagonal",
    "find_common_type",
    "finfo",
    "fix",
    "flatiter",
    "flatnonzero",
    "flexible",
    "flip",
    "fliplr",
    "flipud",
    "float128",
    "float16",
    "float32",
    "float64",
    "float_",
    "float_power",
    "floating",
    "floor",
    "floor_divide",
    "fmax",
    "fmin",
    "fmod",
    "format_float_positional",
    "format_float_scientific",
    "format_parser",
    "frexp",
    "from_dlpack",
    "frombuffer",
    "fromfile",
    "fromfunction",
    "fromiter",
    "frompyfunc",
    "fromregex",
    "fromstring",
    "full",
    "full_like",
    "gcd",
    "generic",
    "genfromtxt",
    "geomspace",
    "get_array_wrap",
    "get_include",
    "get_printoptions",
    "getbufsize",
    "geterr",
    "geterrcall",
    "geterrobj",
    "gradient",
    "greater",
    "greater_equal",
    "half",
    "hamming",
    "hanning",
    "heaviside",
    "histogram",
    "histogram2d",
    "histogram_bin_edges",
    "histogramdd",
    "hsplit",
    "hstack",
    "hypot",
    "i0",
    "identity",
    "iinfo",
    "imag",
    "in1d",
    "index_exp",
    "indices",
    "inexact",
    "inf",
    "info",
    "infty",
    "inner",
    "insert",
    "int0",
    "int16",
    "int32",
    "int64",
    "int8",
    "int_",
    "intc",
    "integer",
    "interp",
    "intersect1d",
    "intp",
    "invert",
    "is_busday",
    "isclose",
    "iscomplex",
    "iscomplexobj",
    "isfinite",
    "isfortran",
    "isin",
    "isinf",
    "isnan",
    "isnat",
    "isneginf",
    "isposinf",
    "isreal",
    "isrealobj",
    "isscalar",
    "issctype",
    "issubclass_",
    "issubdtype",
    "issubsctype",
    "iterable",
    "ix_",
    "kaiser",
    "kron",
    "lcm",
    "ldexp",
    "left_shift",
    "less",
    "less_equal",
    "lexsort",
    "linalg",
    "linspace",
    "little_endian",
    "load",
    "loadtxt",
    "log",
    "log10",
    "log1p",
    "log2",
    "logaddexp",
    "logaddexp2",
    "logical_and",
    "logical_not",
    "logical_or",
    "logical_xor",
    "logspace",
    "longcomplex",
    "longdouble",
    "longfloat",
    "longlong",
    "lookfor",
    "ma",
    "mask_indices",
    "mat",
    "math",
    "matmul",
    "matrix",
    "maximum",
    "maximum_sctype",
    "may_share_memory",
    "mean",
    "median",
    "memmap",
    "meshgrid",
    "mgrid",
    "min_scalar_type",
    "minimum",
    "mintypecode",
    "mod",
    "modf",
    "moveaxis",
    "msort",
    "multiply",
    "nan",
    "nan_to_num",
    "nanargmax",
    "nanargmin",
    "nancumprod",
    "nancumsum",
    "nanmax",
    "nanmean",
    "nanmedian",
    "nanmin",
    "nanpercentile",
    "nanprod",
    "nanquantile",
    "nanstd",
    "nansum",
    "nanvar",
    "nbytes",
    "ndarray",
    "ndenumerate",
    "ndim",
    "ndindex",
    "nditer",
    "negative",
    "nested_iters",
    "newaxis",
    "nextafter",
    "nonzero",
    "not_equal",
    "number",
    "obj2sctype",
    "object0",
    "object_",
    "ogrid",
    "ones",
    "ones_like",
    "outer",
    "packbits",
    "pad",
    "partition",
    "percentile",
    "pi",
    "pid",
    "piecewise",
    "place",
    "poly",
    "poly1d",
    "polyadd",
    "polyder",
    "polydiv",
    "polyfit",
    "polyint",
    "polymul",
    "polysub",
    "polyval",
    "positive",
    "power",
    "printoptions",
    "prod",
    "product",
    "promote_types",
    "ptp",
    "put",
    "put_along_axis",
    "putmask",
    "quantile",
    "r_",
    "rad2deg",
    "radians",
    "random",
    "ravel",
    "ravel_multi_index",
    "real",
    "real_if_close",
    "rec",
    "recarray",
    "recfromcsv",
    "recfromtxt",
    "reciprocal",
    "record",
    "remainder",
    "repeat",
    "require",
    "reshape",
    "resize",
    "result_type",
    "right_shift",
    "rint",
    "roll",
    "rollaxis",
    "roots",
    "rot90",
    "round_",
    "row_stack",
    "s_",
    "safe_eval",
    "save",
    "savetxt",
    "savez",
    "savez_compressed",
    "sctype2char",
    "sctypeDict",
    "sctypes",
    "searchsorted",
    "select",
    "set_numeric_ops",
    "set_printoptions",
    "set_string_function",
    "setbufsize",
    "setdiff1d",
    "seterr",
    "seterrcall",
    "seterrobj",
    "setxor1d",
    "shape",
    "shares_memory",
    "short",
    "show_config",
    "sign",
    "signbit",
    "signedinteger",
    "sin",
    "sinc",
    "single",
    "singlecomplex",
    "sinh",
    "size",
    "sometrue",
    "sort",
    "sort_complex",
    "source",
    "spacing",
    "split",
    "sqrt",
    "square",
    "squeeze",
    "stack",
    "std",
    "str0",
    "str_",
    "string_",
    "subtract",
    "sum",
    "swapaxes",
    "take",
    "take_along_axis",
    "tan",
    "tanh",
    "tensordot",
    "tile",
    "timedelta64",
    "trace",
    "tracemalloc_domain",
    "transpose",
    "trapz",
    "tri",
    "tril",
    "tril_indices",
    "tril_indices_from",
    "trim_zeros",
    "triu",
    "triu_indices",
    "triu_indices_from",
    "true_divide",
    "trunc",
    "typecodes",
    "typename",
    "ubyte",
    "ufunc",
    "uint",
    "uint0",
    "uint16",
    "uint32",
    "uint64",
    "uint8",
    "uintc",
    "uintp",
    "ulonglong",
    "unicode_",
    "union1d",
    "unique",
    "unpackbits",
    "unravel_index",
    "unsignedinteger",
    "unwrap",
    "ushort",
    "vander",
    "var",
    "vdot",
    "vectorize",
    "void",
    "void0",
    "vsplit",
    "vstack",
    "where",
    "who",
    "zeros",
    "zeros_like"
]


ALLOW_THREADS = 1
BUFSIZE = 8192
CLIP = 0
ERR_CALL = 3
ERR_DEFAULT = 521
ERR_IGNORE = 0
ERR_LOG = 5
ERR_PRINT = 4
ERR_RAISE = 2
ERR_WARN = 1
FLOATING_POINT_SUPPORT = 1
FPE_DIVIDEBYZERO = 1
FPE_INVALID = 8
FPE_OVERFLOW = 2
FPE_UNDERFLOW = 4
False_: numpy.bool_ # value = False
Inf: float # value = inf
Infinity: float # value = inf
MAXDIMS = 32
MAY_SHARE_BOUNDS = 0
MAY_SHARE_EXACT = -1
NAN: float # value = nan
NINF: float # value = -inf
NZERO = -0.0
NaN: float # value = nan
PINF: float # value = inf
PZERO = 0.0
RAISE = 2
SHIFT_DIVIDEBYZERO = 0
SHIFT_INVALID = 9
SHIFT_OVERFLOW = 3
SHIFT_UNDERFLOW = 6
ScalarType: tuple # value = (<class 'int'>, <class 'float'>, <class 'complex'>, <class 'bool'>, <class 'bytes'>, <class 'str'>, <class 'memoryview'>, <class 'numpy.bool_'>, <class 'numpy.complex64'>, <class 'numpy.complex128'>, <class 'numpy.complex256'>, <class 'numpy.float16'>, <class 'numpy.float32'>, <class 'numpy.float64'>, <class 'numpy.float128'>, <class 'numpy.int8'>, <class 'numpy.int16'>, <class 'numpy.int32'>, <class 'numpy.longlong'>, <class 'numpy.int64'>, <class 'numpy.timedelta64'>, <class 'numpy.datetime64'>, <class 'numpy.object_'>, <class 'numpy.bytes_'>, <class 'numpy.str_'>, <class 'numpy.uint8'>, <class 'numpy.uint16'>, <class 'numpy.uint32'>, <class 'numpy.ulonglong'>, <class 'numpy.uint64'>, <class 'numpy.void'>)
True_: numpy.bool_ # value = True
UFUNC_BUFSIZE_DEFAULT = 8192
UFUNC_PYVALS_NAME = 'UFUNC_PYVALS'
WRAP = 1
_UFUNC_API: typing.Any  # PyCapsule()
__version__ = '1.23.5'
absolute: numpy.ufunc # value = <ufunc 'absolute'>
add: numpy.ufunc # value = <ufunc 'add'>
arccos: numpy.ufunc # value = <ufunc 'arccos'>
arccosh: numpy.ufunc # value = <ufunc 'arccosh'>
arcsin: numpy.ufunc # value = <ufunc 'arcsin'>
arcsinh: numpy.ufunc # value = <ufunc 'arcsinh'>
arctan: numpy.ufunc # value = <ufunc 'arctan'>
arctan2: numpy.ufunc # value = <ufunc 'arctan2'>
arctanh: numpy.ufunc # value = <ufunc 'arctanh'>
bitwise_and: numpy.ufunc # value = <ufunc 'bitwise_and'>
bitwise_not: numpy.ufunc # value = <ufunc 'invert'>
bitwise_or: numpy.ufunc # value = <ufunc 'bitwise_or'>
bitwise_xor: numpy.ufunc # value = <ufunc 'bitwise_xor'>
c_: numpy.lib.index_tricks.CClass
cast: numpy.core.numerictypes._typedict # value = {<class 'numpy.float32'>: <function <lambda>>, <class 'numpy.uint16'>: <function <lambda>>, <class 'numpy.int16'>: <function <lambda>>, <class 'numpy.bool_'>: <function <lambda>>, <class 'numpy.timedelta64'>: <function <lambda>>, <class 'numpy.float16'>: <function <lambda>>, <class 'numpy.uint8'>: <function <lambda>>, <class 'numpy.int8'>: <function <lambda>>, <class 'numpy.object_'>: <function <lambda>>, <class 'numpy.datetime64'>: <function <lambda>>, <class 'numpy.ulonglong'>: <function <lambda>>, <class 'numpy.longlong'>: <function <lambda>>, <class 'numpy.void'>: <function <lambda>>, <class 'numpy.complex256'>: <function <lambda>>, <class 'numpy.float128'>: <function <lambda>>, <class 'numpy.uint64'>: <function <lambda>>, <class 'numpy.int64'>: <function <lambda>>, <class 'numpy.str_'>: <function <lambda>>, <class 'numpy.complex128'>: <function <lambda>>, <class 'numpy.float64'>: <function <lambda>>, <class 'numpy.uint32'>: <function <lambda>>, <class 'numpy.int32'>: <function <lambda>>, <class 'numpy.bytes_'>: <function <lambda>>, <class 'numpy.complex64'>: <function <lambda>>}
cbrt: numpy.ufunc # value = <ufunc 'cbrt'>
ceil: numpy.ufunc # value = <ufunc 'ceil'>
conj: numpy.ufunc # value = <ufunc 'conjugate'>
conjugate: numpy.ufunc # value = <ufunc 'conjugate'>
copysign: numpy.ufunc # value = <ufunc 'copysign'>
cos: numpy.ufunc # value = <ufunc 'cos'>
cosh: numpy.ufunc # value = <ufunc 'cosh'>
deg2rad: numpy.ufunc # value = <ufunc 'deg2rad'>
degrees: numpy.ufunc # value = <ufunc 'degrees'>
divide: numpy.ufunc # value = <ufunc 'divide'>
divmod: numpy.ufunc # value = <ufunc 'divmod'>
e = 2.718281828459045
equal: numpy.ufunc # value = <ufunc 'equal'>
euler_gamma = 0.5772156649015329
exp: numpy.ufunc # value = <ufunc 'exp'>
exp2: numpy.ufunc # value = <ufunc 'exp2'>
expm1: numpy.ufunc # value = <ufunc 'expm1'>
fabs: numpy.ufunc # value = <ufunc 'fabs'>
float_power: numpy.ufunc # value = <ufunc 'float_power'>
floor: numpy.ufunc # value = <ufunc 'floor'>
floor_divide: numpy.ufunc # value = <ufunc 'floor_divide'>
fmax: numpy.ufunc # value = <ufunc 'fmax'>
fmin: numpy.ufunc # value = <ufunc 'fmin'>
fmod: numpy.ufunc # value = <ufunc 'fmod'>
frexp: numpy.ufunc # value = <ufunc 'frexp'>
gcd: numpy.ufunc # value = <ufunc 'gcd'>
greater: numpy.ufunc # value = <ufunc 'greater'>
greater_equal: numpy.ufunc # value = <ufunc 'greater_equal'>
heaviside: numpy.ufunc # value = <ufunc 'heaviside'>
hypot: numpy.ufunc # value = <ufunc 'hypot'>
index_exp: numpy.lib.index_tricks.IndexExpression
inf: float # value = inf
infty: float # value = inf
invert: numpy.ufunc # value = <ufunc 'invert'>
isfinite: numpy.ufunc # value = <ufunc 'isfinite'>
isinf: numpy.ufunc # value = <ufunc 'isinf'>
isnan: numpy.ufunc # value = <ufunc 'isnan'>
isnat: numpy.ufunc # value = <ufunc 'isnat'>
lcm: numpy.ufunc # value = <ufunc 'lcm'>
ldexp: numpy.ufunc # value = <ufunc 'ldexp'>
left_shift: numpy.ufunc # value = <ufunc 'left_shift'>
less: numpy.ufunc # value = <ufunc 'less'>
less_equal: numpy.ufunc # value = <ufunc 'less_equal'>
little_endian = True
log: numpy.ufunc # value = <ufunc 'log'>
log10: numpy.ufunc # value = <ufunc 'log10'>
log1p: numpy.ufunc # value = <ufunc 'log1p'>
log2: numpy.ufunc # value = <ufunc 'log2'>
logaddexp: numpy.ufunc # value = <ufunc 'logaddexp'>
logaddexp2: numpy.ufunc # value = <ufunc 'logaddexp2'>
logical_and: numpy.ufunc # value = <ufunc 'logical_and'>
logical_not: numpy.ufunc # value = <ufunc 'logical_not'>
logical_or: numpy.ufunc # value = <ufunc 'logical_or'>
logical_xor: numpy.ufunc # value = <ufunc 'logical_xor'>
matmul: numpy.ufunc # value = <ufunc 'matmul'>
maximum: numpy.ufunc # value = <ufunc 'maximum'>
mgrid: numpy.lib.index_tricks.MGridClass
minimum: numpy.ufunc # value = <ufunc 'minimum'>
mod: numpy.ufunc # value = <ufunc 'remainder'>
modf: numpy.ufunc # value = <ufunc 'modf'>
multiply: numpy.ufunc # value = <ufunc 'multiply'>
nan: float # value = nan
nbytes: numpy.core.numerictypes._typedict # value = {<class 'numpy.bool_'>: 1, <class 'numpy.int8'>: 1, <class 'numpy.uint8'>: 1, <class 'numpy.int16'>: 2, <class 'numpy.uint16'>: 2, <class 'numpy.int32'>: 4, <class 'numpy.uint32'>: 4, <class 'numpy.int64'>: 8, <class 'numpy.uint64'>: 8, <class 'numpy.longlong'>: 8, <class 'numpy.ulonglong'>: 8, <class 'numpy.float16'>: 2, <class 'numpy.float32'>: 4, <class 'numpy.float64'>: 8, <class 'numpy.float128'>: 16, <class 'numpy.complex64'>: 8, <class 'numpy.complex128'>: 16, <class 'numpy.complex256'>: 32, <class 'numpy.object_'>: 8, <class 'numpy.bytes_'>: 0, <class 'numpy.str_'>: 0, <class 'numpy.void'>: 0, <class 'numpy.datetime64'>: 8, <class 'numpy.timedelta64'>: 8}
negative: numpy.ufunc # value = <ufunc 'negative'>
newaxis = None
nextafter: numpy.ufunc # value = <ufunc 'nextafter'>
not_equal: numpy.ufunc # value = <ufunc 'not_equal'>
ogrid: numpy.lib.index_tricks.OGridClass
pi = 3.141592653589793
positive: numpy.ufunc # value = <ufunc 'positive'>
power: numpy.ufunc # value = <ufunc 'power'>
r_: numpy.lib.index_tricks.RClass
rad2deg: numpy.ufunc # value = <ufunc 'rad2deg'>
radians: numpy.ufunc # value = <ufunc 'radians'>
reciprocal: numpy.ufunc # value = <ufunc 'reciprocal'>
remainder: numpy.ufunc # value = <ufunc 'remainder'>
right_shift: numpy.ufunc # value = <ufunc 'right_shift'>
rint: numpy.ufunc # value = <ufunc 'rint'>
s_: numpy.lib.index_tricks.IndexExpression
sctypeDict: dict # value = {'?': <class 'numpy.bool_'>, 0: <class 'numpy.bool_'>, 'byte': <class 'numpy.int8'>, 'b': <class 'numpy.int8'>, 1: <class 'numpy.int8'>, 'ubyte': <class 'numpy.uint8'>, 'B': <class 'numpy.uint8'>, 2: <class 'numpy.uint8'>, 'short': <class 'numpy.int16'>, 'h': <class 'numpy.int16'>, 3: <class 'numpy.int16'>, 'ushort': <class 'numpy.uint16'>, 'H': <class 'numpy.uint16'>, 4: <class 'numpy.uint16'>, 'i': <class 'numpy.int32'>, 5: <class 'numpy.int32'>, 'uint': <class 'numpy.uint64'>, 'I': <class 'numpy.uint32'>, 6: <class 'numpy.uint32'>, 'intp': <class 'numpy.int64'>, 'p': <class 'numpy.int64'>, 7: <class 'numpy.int64'>, 'uintp': <class 'numpy.uint64'>, 'P': <class 'numpy.uint64'>, 8: <class 'numpy.uint64'>, 'long': <class 'numpy.int64'>, 'l': <class 'numpy.int64'>, 'ulong': <class 'numpy.uint64'>, 'L': <class 'numpy.uint64'>, 'longlong': <class 'numpy.longlong'>, 'q': <class 'numpy.longlong'>, 9: <class 'numpy.longlong'>, 'ulonglong': <class 'numpy.ulonglong'>, 'Q': <class 'numpy.ulonglong'>, 10: <class 'numpy.ulonglong'>, 'half': <class 'numpy.float16'>, 'e': <class 'numpy.float16'>, 23: <class 'numpy.float16'>, 'f': <class 'numpy.float32'>, 11: <class 'numpy.float32'>, 'double': <class 'numpy.float64'>, 'd': <class 'numpy.float64'>, 12: <class 'numpy.float64'>, 'longdouble': <class 'numpy.float128'>, 'g': <class 'numpy.float128'>, 13: <class 'numpy.float128'>, 'cfloat': <class 'numpy.complex128'>, 'F': <class 'numpy.complex64'>, 14: <class 'numpy.complex64'>, 'cdouble': <class 'numpy.complex128'>, 'D': <class 'numpy.complex128'>, 15: <class 'numpy.complex128'>, 'clongdouble': <class 'numpy.complex256'>, 'G': <class 'numpy.complex256'>, 16: <class 'numpy.complex256'>, 'O': <class 'numpy.object_'>, 17: <class 'numpy.object_'>, 'S': <class 'numpy.bytes_'>, 18: <class 'numpy.bytes_'>, 'unicode': <class 'numpy.str_'>, 'U': <class 'numpy.str_'>, 19: <class 'numpy.str_'>, 'void': <class 'numpy.void'>, 'V': <class 'numpy.void'>, 20: <class 'numpy.void'>, 'M': <class 'numpy.datetime64'>, 21: <class 'numpy.datetime64'>, 'm': <class 'numpy.timedelta64'>, 22: <class 'numpy.timedelta64'>, 'bool8': <class 'numpy.bool_'>, 'b1': <class 'numpy.bool_'>, 'int64': <class 'numpy.int64'>, 'i8': <class 'numpy.int64'>, 'uint64': <class 'numpy.uint64'>, 'u8': <class 'numpy.uint64'>, 'float16': <class 'numpy.float16'>, 'f2': <class 'numpy.float16'>, 'float32': <class 'numpy.float32'>, 'f4': <class 'numpy.float32'>, 'float64': <class 'numpy.float64'>, 'f8': <class 'numpy.float64'>, 'float128': <class 'numpy.float128'>, 'f16': <class 'numpy.float128'>, 'complex64': <class 'numpy.complex64'>, 'c8': <class 'numpy.complex64'>, 'complex128': <class 'numpy.complex128'>, 'c16': <class 'numpy.complex128'>, 'complex256': <class 'numpy.complex256'>, 'c32': <class 'numpy.complex256'>, 'object0': <class 'numpy.object_'>, 'bytes0': <class 'numpy.bytes_'>, 'str0': <class 'numpy.str_'>, 'void0': <class 'numpy.void'>, 'datetime64': <class 'numpy.datetime64'>, 'M8': <class 'numpy.datetime64'>, 'timedelta64': <class 'numpy.timedelta64'>, 'm8': <class 'numpy.timedelta64'>, 'int32': <class 'numpy.int32'>, 'i4': <class 'numpy.int32'>, 'uint32': <class 'numpy.uint32'>, 'u4': <class 'numpy.uint32'>, 'int16': <class 'numpy.int16'>, 'i2': <class 'numpy.int16'>, 'uint16': <class 'numpy.uint16'>, 'u2': <class 'numpy.uint16'>, 'int8': <class 'numpy.int8'>, 'i1': <class 'numpy.int8'>, 'uint8': <class 'numpy.uint8'>, 'u1': <class 'numpy.uint8'>, 'complex_': <class 'numpy.complex128'>, 'int0': <class 'numpy.int64'>, 'uint0': <class 'numpy.uint64'>, 'single': <class 'numpy.float32'>, 'csingle': <class 'numpy.complex64'>, 'singlecomplex': <class 'numpy.complex64'>, 'float_': <class 'numpy.float64'>, 'intc': <class 'numpy.int32'>, 'uintc': <class 'numpy.uint32'>, 'int_': <class 'numpy.int64'>, 'longfloat': <class 'numpy.float128'>, 'clongfloat': <class 'numpy.complex256'>, 'longcomplex': <class 'numpy.complex256'>, 'bool_': <class 'numpy.bool_'>, 'bytes_': <class 'numpy.bytes_'>, 'string_': <class 'numpy.bytes_'>, 'str_': <class 'numpy.str_'>, 'unicode_': <class 'numpy.str_'>, 'object_': <class 'numpy.object_'>, 'int': <class 'numpy.int64'>, 'float': <class 'numpy.float64'>, 'complex': <class 'numpy.complex128'>, 'bool': <class 'numpy.bool_'>, 'object': <class 'numpy.object_'>, 'str': <class 'numpy.str_'>, 'bytes': <class 'numpy.bytes_'>, 'a': <class 'numpy.bytes_'>}
sctypes: dict # value = {'int': [<class 'numpy.int8'>, <class 'numpy.int16'>, <class 'numpy.int32'>, <class 'numpy.int64'>], 'uint': [<class 'numpy.uint8'>, <class 'numpy.uint16'>, <class 'numpy.uint32'>, <class 'numpy.uint64'>], 'float': [<class 'numpy.float16'>, <class 'numpy.float32'>, <class 'numpy.float64'>, <class 'numpy.float128'>], 'complex': [<class 'numpy.complex64'>, <class 'numpy.complex128'>, <class 'numpy.complex256'>], 'others': [<class 'bool'>, <class 'object'>, <class 'bytes'>, <class 'str'>, <class 'numpy.void'>]}
sign: numpy.ufunc # value = <ufunc 'sign'>
signbit: numpy.ufunc # value = <ufunc 'signbit'>
sin: numpy.ufunc # value = <ufunc 'sin'>
sinh: numpy.ufunc # value = <ufunc 'sinh'>
spacing: numpy.ufunc # value = <ufunc 'spacing'>
sqrt: numpy.ufunc # value = <ufunc 'sqrt'>
square: numpy.ufunc # value = <ufunc 'square'>
subtract: numpy.ufunc # value = <ufunc 'subtract'>
tan: numpy.ufunc # value = <ufunc 'tan'>
tanh: numpy.ufunc # value = <ufunc 'tanh'>
tracemalloc_domain = 389047
true_divide: numpy.ufunc # value = <ufunc 'divide'>
trunc: numpy.ufunc # value = <ufunc 'trunc'>
typecodes = {'Character': 'c', 'Integer': 'bhilqp', 'UnsignedInteger': 'BHILQP', 'Float': 'efdg', 'Complex': 'FDG', 'AllInteger': 'bBhHiIlLqQpP', 'AllFloat': 'efdgFDG', 'Datetime': 'Mm', 'All': '?bhilqpBHILQPefdgFDGSUVOMm'}
bool8 = numpy.bool_
byte = numpy.int8
bytes0 = numpy.bytes_
cdouble = numpy.complex128
cfloat = numpy.complex128
clongdouble = numpy.complex256
clongfloat = numpy.complex256
complex_ = numpy.complex128
csingle = numpy.complex64
double = numpy.float64
float_ = numpy.float64
half = numpy.float16
int0 = numpy.int64
int_ = numpy.int64
intc = numpy.int32
intp = numpy.int64
longcomplex = numpy.complex256
longdouble = numpy.float128
longfloat = numpy.float128
mat = numpy.asmatrix
object0 = numpy.object_
row_stack = numpy.vstack
short = numpy.int16
show_config = numpy.__config__.show
single = numpy.float32
singlecomplex = numpy.complex64
str0 = numpy.str_
string_ = numpy.bytes_
ubyte = numpy.uint8
uint = numpy.uint64
uint0 = numpy.uint64
uintc = numpy.uint32
uintp = numpy.uint64
unicode_ = numpy.str_
ushort = numpy.uint16
void0 = numpy.void
