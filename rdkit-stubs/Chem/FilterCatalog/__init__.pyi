from __future__ import annotations
import rdkit.Chem.FilterCatalog
import typing
from rdkit.Chem.rdfiltercatalog import ExclusionList
from rdkit.Chem.rdfiltercatalog import FilterCatalog
from rdkit.Chem.rdfiltercatalog import FilterCatalogEntry
from rdkit.Chem.rdfiltercatalog import FilterCatalogEntryList
from rdkit.Chem.rdfiltercatalog import FilterCatalogListOfEntryList
from rdkit.Chem.rdfiltercatalog import FilterCatalogParams
from rdkit.Chem.rdfiltercatalog import FilterHierarchyMatcher
from rdkit.Chem.rdfiltercatalog import FilterMatch
from rdkit.Chem.rdfiltercatalog import FilterMatcherBase
from rdkit.Chem.rdfiltercatalog import IntPair
from rdkit.Chem.rdfiltercatalog import MatchTypeVect
from rdkit.Chem.rdfiltercatalog import MolList
from rdkit.Chem.rdfiltercatalog import PythonFilterMatcher
from rdkit.Chem.rdfiltercatalog import SmartsMatcher
from rdkit.Chem.rdfiltercatalog import VectFilterMatch
import Boost.Python
import rdkit.Chem
import rdkit.Chem.rdfiltercatalog
import rdkit.Chem.rdfiltercatalog.FilterMatchOps
import sys

__all__ = [
    "Chem",
    "ExclusionList",
    "FilterCatalog",
    "FilterCatalogCanSerialize",
    "FilterCatalogEntry",
    "FilterCatalogEntryList",
    "FilterCatalogListOfEntryList",
    "FilterCatalogParams",
    "FilterHierarchyMatcher",
    "FilterMatch",
    "FilterMatchOps",
    "FilterMatcher",
    "FilterMatcherBase",
    "GetFlattenedFunctionalGroupHierarchy",
    "GetFunctionalGroupHierarchy",
    "IntPair",
    "MatchTypeVect",
    "MolList",
    "PythonFilterMatcher",
    "RunFilterCatalog",
    "SmartsMatcher",
    "VectFilterMatch",
    "sys"
]


class FilterMatcher(rdkit.Chem.rdfiltercatalog.PythonFilterMatcher, rdkit.Chem.rdfiltercatalog.FilterMatcherBase, Boost.Python.instance):
    """
    FilterMatcher - This class allows creation of Python based
        filters.  Subclass this class to create a Filter useable
        in a FilterCatalogEntry

        Simple Example:

        from rdkit.Chem import rdMolDescriptors
        class MWFilter(FilterMatcher):
          def __init__(self, minMw, maxMw):
              FilterMatcher.__init__(self, "MW violation")
              self.minMw = minMw
              self.maxMw = maxMw

          def IsValid(self):
             return True

          def HasMatch(self, mol):
             mw = rdMolDescriptors.CalcExactMolWt(mol)
             return not self.minMw <= mw <= self.maxMw
        
    """
    pass
def FilterCatalogCanSerialize() -> bool:
    """
    FilterCatalogCanSerialize() -> bool
        Returns True if the FilterCatalog is serializable (requires boost serialization

        C++ signature :
            bool FilterCatalogCanSerialize()
    """
def GetFlattenedFunctionalGroupHierarchy(normalized: bool = False) -> dict:
    """
    GetFlattenedFunctionalGroupHierarchy(normalized: bool = False) -> dict
        Returns the flattened functional group hierarchy as a dictionary  of name:ROMOL_SPTR substructure items

        C++ signature :
            boost::python::dict GetFlattenedFunctionalGroupHierarchy([ bool=False])
    """
def GetFunctionalGroupHierarchy() -> FilterCatalog:
    """
    GetFunctionalGroupHierarchy() -> FilterCatalog
        Returns the functional group hierarchy filter catalog

        C++ signature :
            RDKit::FilterCatalog GetFunctionalGroupHierarchy()
    """
def RunFilterCatalog( filterCatalog: FilterCatalog, smiles: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, numThreads: int = 1) -> FilterCatalogListOfEntryList:
    """
    RunFilterCatalog( filterCatalog: FilterCatalog, smiles: _vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE, numThreads: int = 1) -> FilterCatalogListOfEntryList
        Run the filter catalog on the input list of smiles strings.
        Use numThreads=0 to use all available processors. Returns a vector of vectors.  For each input smiles, a vector of FilterCatalogEntry objects are returned for each matched filter.  If a molecule matches no filter, the vector will be empty. If a smiles string can't be parsed, a 'Bad smiles' entry is returned.

        C++ signature :
            std::vector<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > >, std::allocator<std::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const> > > > > RunFilterCatalog(RDKit::FilterCatalog,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > [,int=1])
    """
