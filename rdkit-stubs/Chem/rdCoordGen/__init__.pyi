"""Module containing interface to the CoordGen library."""
from __future__ import annotations
import rdkit.Chem.rdCoordGen
import typing
import Boost.Python

__all__ = [
    "AddCoords",
    "CoordGenParams",
    "SetDefaultTemplateFileDir"
]


class CoordGenParams(Boost.Python.instance):
    """
    Parameters controlling coordinate generation
    """
    @staticmethod
    def SetCoordMap( arg1: CoordGenParams, arg2: dict) -> None: 
        """
        SetCoordMap( arg1: CoordGenParams, arg2: dict) -> None
            expects a dictionary of Point2D objects with template coordinates

            C++ signature :
                void SetCoordMap(RDKit::CoordGen::CoordGenParams*,boost::python::dict {lvalue})
        """
    @staticmethod
    def SetTemplateMol( arg1: CoordGenParams, arg2: Mol) -> None: 
        """
        SetTemplateMol( arg1: CoordGenParams, arg2: Mol) -> None
            sets a molecule to be used as the template

            C++ signature :
                void SetTemplateMol(RDKit::CoordGen::CoordGenParams*,RDKit::ROMol const*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def coordgenScaling(self) -> None:
        """
        scaling factor for a single bond

        :type: None
        """
    @property
    def dbg_useConstrained(self) -> None:
        """
        for debugging use

        :type: None
        """
    @property
    def dbg_useFixed(self) -> None:
        """
        for debugging use

        :type: None
        """
    @property
    def minimizerPrecision(self) -> None:
        """
        controls sketcher precision

        :type: None
        """
    @property
    def sketcherBestPrecision(self) -> None:
        """
        highest quality (and slowest) precision setting

        :type: None
        """
    @property
    def sketcherCoarsePrecision(self) -> None:
        """
        "coarse" (fastest) precision setting, produces good-quality coordinates most of the time, this is the default setting for the RDKit

        :type: None
        """
    @property
    def sketcherQuickPrecision(self) -> None:
        """
        faster precision setting

        :type: None
        """
    @property
    def sketcherStandardPrecision(self) -> None:
        """
        standard quality precision setting, the default for the coordgen project

        :type: None
        """
    @property
    def templateFileDir(self) -> None:
        """
        directory containing the templates.mae file

        :type: None
        """
    @property
    def treatNonterminalBondsToMetalAsZOBs(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 144
    pass
def AddCoords( mol: Mol, params: object = None) -> None:
    """
    AddCoords( mol: Mol, params: object = None) -> None
        Add 2D coordinates.
        ARGUMENTS:
           - mol: molecule to modify
           - params: (optional) parameters controlling the coordinate generation
        
        

        C++ signature :
            void AddCoords(RDKit::ROMol {lvalue} [,boost::python::api::object {lvalue}=None])
    """
def SetDefaultTemplateFileDir( arg1: str) -> None:
    """
    SetDefaultTemplateFileDir( arg1: str) -> None

        C++ signature :
            void SetDefaultTemplateFileDir(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
