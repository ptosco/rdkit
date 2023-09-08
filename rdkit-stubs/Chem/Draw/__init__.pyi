from __future__ import annotations
import rdkit.Chem.Draw
import typing
from _io import BytesIO
from rdkit.Chem.Draw.rdMolDraw2D import ContourParams
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from rdkit.Chem.Draw.rdMolDraw2D import IntStringMap
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2D
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DCairo
from rdkit.Chem.Draw.rdMolDraw2D import MolDraw2DSVG
from rdkit.Chem.Draw.rdMolDraw2D import MolDrawOptions
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem.Draw.rdMolDraw2D import MultiColourHighlightStyle
from rdkit.Chem.Draw.rdMolDraw2D import map_indexing_suite_IntStringMap_entry
import _collections
import numpy
import os
import rdkit.Chem
import rdkit.Chem.Draw.rdMolDraw2D
import rdkit.Chem.rdDepictor
import rdkit.RDConfig
import rdkit.rdBase
import warnings
_Shape = typing.Tuple[int, ...]

__all__ = [
    "BytesIO",
    "Chem",
    "CircleAndLine",
    "ContourAndDrawGaussians",
    "ContourAndDrawGrid",
    "ContourParams",
    "DebugDraw",
    "DrawMoleculeACS1996",
    "DrawMorganBit",
    "DrawMorganBits",
    "DrawMorganEnv",
    "DrawMorganEnvs",
    "DrawRDKitBit",
    "DrawRDKitBits",
    "DrawRDKitEnv",
    "DrawRDKitEnvs",
    "DrawingOptions",
    "FingerprintEnv",
    "IntStringMap",
    "Lasso",
    "MeanBondLength",
    "MolDraw2D",
    "MolDraw2DCairo",
    "MolDraw2DSVG",
    "MolDrawOptions",
    "MolDrawing",
    "MolToACS1996SVG",
    "MolToFile",
    "MolToImage",
    "MolToImageFile",
    "MolToMPL",
    "MolToQPixmap",
    "MolToSVG",
    "MolsMatrixToGridImage",
    "MolsToGridImage",
    "MolsToImage",
    "MultiColourHighlightStyle",
    "PrepareAndDrawMolecule",
    "PrepareMolForDrawing",
    "RDConfig",
    "ReactionToImage",
    "SetACS1996Mode",
    "SetComicMode",
    "SetDarkMode",
    "SetMonochromeMode",
    "ShowMol",
    "UpdateDrawerParamsFromJSON",
    "UpdateMolDrawOptionsFromJSON",
    "calcAtomGaussians",
    "find_spec",
    "map_indexing_suite_IntStringMap_entry",
    "namedtuple",
    "numpy",
    "os",
    "rdBase",
    "rdDepictor",
    "rdMolDraw2D",
    "shouldKekulize",
    "warnings"
]


class FingerprintEnv(tuple):
    """
    FingerprintEnv(submol, highlightAtoms, atomColors, highlightBonds, bondColors, highlightRadii)
    """
    __slots__ = ()
    _field_defaults = {}
    _fields = ('submol', 'highlightAtoms', 'atomColors', 'highlightBonds', 'bondColors', 'highlightRadii')
    _fields_defaults = {}
    atomColors: _collections._tuplegetter
    bondColors: _collections._tuplegetter
    highlightAtoms: _collections._tuplegetter
    highlightBonds: _collections._tuplegetter
    highlightRadii: _collections._tuplegetter
    submol: _collections._tuplegetter
    pass
def ContourAndDrawGaussians( drawer: MolDraw2D, locs: AtomPairsParameters, heights: AtomPairsParameters, widths: AtomPairsParameters, nContours: int = 10, levels: AtomPairsParameters = None, params: ContourParams = ContourParams(), mol: AtomPairsParameters = None) -> None:
    """
    ContourAndDrawGaussians( drawer: MolDraw2D, locs: AtomPairsParameters, heights: AtomPairsParameters, widths: AtomPairsParameters, nContours: int = 10, levels: AtomPairsParameters = None, params: ContourParams = ContourParams(), mol: AtomPairsParameters = None) -> None
        Generates and draws contours for a set of gaussians
        
          - drawer: the MolDraw2D object to use
          - locs: locations of the gaussians
          - heights: the heights (or weights) of the gaussians
          - widths: the standard deviations of the gaussians
          - nContours: the number of contours to draw
          - levels: the contours to use
          - ps: additional parameters controlling the contouring.
          - mol: molecule used to help set scale.
        
          The values are calculated on a grid with spacing params.gridResolution.
          If params.setScale  is set, the grid size will be calculated based on the
          locations of the gaussians and params.extraGridPadding. Otherwise the current
          size of the viewport will be used.
        
          If the levels argument is empty, the contour levels will be determined
          automatically from the max and min values on the grid and levels will
          be updated to include the contour levels.
        
          If params.fillGrid is set, the data on the grid will also be drawn using
          the color scheme in params.colourMap
        
          If mol is not 0, uses the molecule to help set the scale, assuming that
          it will be drawn over the plot, so needs to fit on it.
        */

        C++ signature :
            void ContourAndDrawGaussians(RDKit::MolDraw2D {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object [,unsigned int=10 [,boost::python::api::object=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7f2c6425a640> [,boost::python::api::object=None]]]])
    """
def ContourAndDrawGrid( drawer: MolDraw2D, data: AtomPairsParameters, xcoords: AtomPairsParameters, ycoords: AtomPairsParameters, nContours: int = 10, levels: AtomPairsParameters = None, params: ContourParams = ContourParams(), mol: AtomPairsParameters = None) -> None:
    """
    ContourAndDrawGrid( drawer: MolDraw2D, data: AtomPairsParameters, xcoords: AtomPairsParameters, ycoords: AtomPairsParameters, nContours: int = 10, levels: AtomPairsParameters = None, params: ContourParams = ContourParams(), mol: AtomPairsParameters = None) -> None
        Generates and draws contours for data on a grid
        
          - drawer: the MolDraw2D object to use
          - data: numpy array with the data to be contoured
          - xcoords: the x coordinates of the grid
          - ycoords: the y coordinates of the grid
          - nContours: the number of contours to draw
          - levels: the contours to use
          - ps: additional parameters controlling the contouring
          - mol: molecule used to help set scale.
        
          The values are calculated on a grid with spacing params.gridResolution.
          If params.setScale  is set, the grid size will be calculated based on the
          locations of the gaussians and params.extraGridPadding. Otherwise the current
          size of the viewport will be used.
        
          If the levels argument is empty, the contour levels will be determined
          automatically from the max and min values on the grid and levels will
          be updated to include the contour levels.
        
          If params.fillGrid is set, the data on the grid will also be drawn using
          the color scheme in params.colourMap
        
          If mol is not 0, uses the molecule to help set the scale, assuming that
          it will be drawn over the plot, so needs to fit on it.
        */

        C++ signature :
            void ContourAndDrawGrid(RDKit::MolDraw2D {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue} [,unsigned int=10 [,boost::python::api::object {lvalue}=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7f2c6425a700> [,boost::python::api::object=None]]]])
    """
def DrawMoleculeACS1996( drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1) -> None:
    """
    DrawMoleculeACS1996( drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1) -> None
        Draws molecule in ACS 1996 mode.

        C++ signature :
            void DrawMoleculeACS1996(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1]]]]]]])
    """
def MeanBondLength( mol: Mol, confId: int = -1) -> float:
    """
    MeanBondLength( mol: Mol, confId: int = -1) -> float
        Calculate the mean bond length for the molecule.

        C++ signature :
            double MeanBondLength(RDKit::ROMol [,int=-1])
    """
def MolToACS1996SVG( mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1) -> str:
    """
    MolToACS1996SVG( mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1) -> str
        Returns ACS 1996 mode svg for a molecule

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToACS1996SVG(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1]]]]]]])
    """
def MolToSVG( mol: Mol, width: int = 300, height: int = 300, highlightAtoms: AtomPairsParameters = None, kekulize: bool = True, lineWidthMult: int = 1, fontSize: bool = 12, includeAtomCircles: int = True) -> str:
    """
    MolToSVG( mol: Mol, width: int = 300, height: int = 300, highlightAtoms: AtomPairsParameters = None, kekulize: bool = True, lineWidthMult: int = 1, fontSize: bool = 12, includeAtomCircles: int = True) -> str
        Returns svg for a molecule

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToSVG(RDKit::ROMol [,unsigned int=300 [,unsigned int=300 [,boost::python::api::object=None [,bool=True [,unsigned int=1 [,bool=12 [,int=True]]]]]]])
    """
def PrepareAndDrawMolecule( drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, kekulize: bool = True) -> None:
    """
    PrepareAndDrawMolecule( drawer: MolDraw2D, mol: Mol, legend: str = '', highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, kekulize: bool = True) -> None
        Preps a molecule for drawing and actually draws it
        

        C++ signature :
            void PrepareAndDrawMolecule(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,bool=True]]]]]]]])
    """
def PrepareMolForDrawing( mol: Mol, kekulize: bool = True, addChiralHs: bool = True, wedgeBonds: bool = True, forceCoords: bool = False, wavyBonds: bool = False) -> Mol:
    """
    PrepareMolForDrawing( mol: Mol, kekulize: bool = True, addChiralHs: bool = True, wedgeBonds: bool = True, forceCoords: bool = False, wavyBonds: bool = False) -> Mol
        Does some cleanup operations on the molecule to prepare it to draw nicely.
        The operations include: kekulization, addition of chiral Hs (so that we can draw
        wedges to them), wedging of bonds at chiral centers, and generation of a 2D
        conformation if the molecule does not already have a conformation
        
        Returns a modified copy of the molecule.
        

        C++ signature :
            RDKit::ROMol* PrepareMolForDrawing(RDKit::ROMol const* [,bool=True [,bool=True [,bool=True [,bool=False [,bool=False]]]]])
    """
def SetACS1996Mode( drawOptions: MolDrawOptions, meanBondLength: float) -> None:
    """
    SetACS1996Mode( drawOptions: MolDrawOptions, meanBondLength: float) -> None
        Set the draw options to produce something as close as possible to
        the ACS 1996 guidelines as described at
        https://en.wikipedia.org/wiki/Wikipedia:Manual_of_Style/Chemistry/Structure_drawing
        
         - MolDrawOptions opt - the options what will be changed
         - float meanBondLength - mean bond length of the molecule
        
         Works best if the MolDraw2D object is created with width and height -1 (a
         flexiCanvas).
         The mean bond length may be calculated with MeanBondLength.
         It is used to calculate the offset for the lines in multiple bonds.
        
         Options changed are:
           bondLineWidth = 0.6
           scaleBondWidth = false
           scalingFactor = 14.4 / meanBondLen
           multipleBondOffset = 0.18
           highlightBondWidthMultiplier = 32
           setMonochromeMode - black and white
           fixedFontSize = 10
           additionalAtomLabelPadding = 0.066
           fontFile - if it isn't set already, then if RDBASE is set and the file
                      exists, uses $RDBASE/Data/Fonts/FreeSans.ttf.  Otherwise uses
                      BuiltinRobotoRegular.
         */
        

        C++ signature :
            void SetACS1996Mode(RDKit::MolDrawOptions {lvalue},double)
    """
@typing.overload
def SetDarkMode( arg1: MolDrawOptions) -> None:
    """
    SetDarkMode( arg1: MolDrawOptions) -> None
        set dark mode for a MolDrawOptions object

        C++ signature :
            void SetDarkMode(RDKit::MolDrawOptions {lvalue})

        C++ signature :
            void SetDarkMode(RDKit::MolDraw2D {lvalue})
    """
@typing.overload
def SetDarkMode( arg1: MolDraw2D) -> None:
    pass
@typing.overload
def SetMonochromeMode( options: MolDrawOptions, fgColour: tuple, bgColour: tuple) -> None:
    """
    SetMonochromeMode( options: MolDrawOptions, fgColour: tuple, bgColour: tuple) -> None
        set monochrome mode for a MolDrawOptions object

        C++ signature :
            void SetMonochromeMode(RDKit::MolDrawOptions {lvalue},boost::python::tuple,boost::python::tuple)

        C++ signature :
            void SetMonochromeMode(RDKit::MolDraw2D {lvalue},boost::python::tuple,boost::python::tuple)
    """
@typing.overload
def SetMonochromeMode( drawer: MolDraw2D, fgColour: tuple, bgColour: tuple) -> None:
    pass
def UpdateDrawerParamsFromJSON( drawer: MolDraw2D, json: str) -> None:
    """
    UpdateDrawerParamsFromJSON( drawer: MolDraw2D, json: str) -> None

        C++ signature :
            void UpdateDrawerParamsFromJSON(RDKit::MolDraw2D {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def UpdateMolDrawOptionsFromJSON( opts: MolDrawOptions, json: str) -> None:
    """
    UpdateMolDrawOptionsFromJSON( opts: MolDrawOptions, json: str) -> None

        C++ signature :
            void UpdateMolDrawOptionsFromJSON(RDKit::MolDrawOptions {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
CircleAndLine = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine
Lasso = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso
