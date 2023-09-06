from __future__ import annotations
import rdkit.Chem.rdfiltercatalog.FilterMatchOps
import typing
import Boost.Python
import rdkit.Chem.rdfiltercatalog

__all__ = [
    "And",
    "Not",
    "Or"
]


class And(rdkit.Chem.rdfiltercatalog.FilterMatcherBase, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object, arg2: FilterMatcherBase, arg3: FilterMatcherBase) -> None: 
        """
        __init__( arg1: object, arg2: FilterMatcherBase, arg3: FilterMatcherBase) -> None

            C++ signature :
                void __init__(_object*,RDKit::FilterMatcherBase {lvalue},RDKit::FilterMatcherBase {lvalue})
        """
    __instance_size__ = 32
    pass
class Not(rdkit.Chem.rdfiltercatalog.FilterMatcherBase, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object, arg2: FilterMatcherBase) -> None: 
        """
        __init__( arg1: object, arg2: FilterMatcherBase) -> None

            C++ signature :
                void __init__(_object*,RDKit::FilterMatcherBase {lvalue})
        """
    __instance_size__ = 32
    pass
class Or(rdkit.Chem.rdfiltercatalog.FilterMatcherBase, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object, arg2: FilterMatcherBase, arg3: FilterMatcherBase) -> None: 
        """
        __init__( arg1: object, arg2: FilterMatcherBase, arg3: FilterMatcherBase) -> None

            C++ signature :
                void __init__(_object*,RDKit::FilterMatcherBase {lvalue},RDKit::FilterMatcherBase {lvalue})
        """
    __instance_size__ = 32
    pass
