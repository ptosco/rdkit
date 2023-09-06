from __future__ import annotations
import rdkit.Chem.MolCatalog
import typing
from rdkit.Chem.rdMolCatalog import MolCatalog
from rdkit.Chem.rdMolCatalog import MolCatalogEntry

__all__ = [
    "CreateMolCatalog",
    "MolCatalog",
    "MolCatalogEntry"
]


def CreateMolCatalog() -> MolCatalog:
    """
    CreateMolCatalog() -> MolCatalog

        C++ signature :
            RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>* CreateMolCatalog()
    """
