"""Module containing a C++ implementation of 2D molecule drawing"""
from __future__ import annotations
import rdkit.Chem.Draw.rdMolDraw2D
import typing
import Boost.Python

__all__ = [
    "CircleAndLine",
    "ContourAndDrawGaussians",
    "ContourAndDrawGrid",
    "ContourParams",
    "DrawMoleculeACS1996",
    "IntStringMap",
    "Lasso",
    "MeanBondLength",
    "MolDraw2D",
    "MolDraw2DCairo",
    "MolDraw2DSVG",
    "MolDrawOptions",
    "MolToACS1996SVG",
    "MolToSVG",
    "MultiColourHighlightStyle",
    "PrepareAndDrawMolecule",
    "PrepareMolForDrawing",
    "SetACS1996Mode",
    "SetDarkMode",
    "SetMonochromeMode",
    "UpdateDrawerParamsFromJSON",
    "UpdateMolDrawOptionsFromJSON",
    "map_indexing_suite_IntStringMap_entry"
]


class ContourParams(Boost.Python.instance):
    """
    Parameters for drawing contours
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def setColourMap(self, colours: AtomPairsParameters) -> None: 
        """
        setColourMap( self: ContourParams, colours: AtomPairsParameters) -> None

            C++ signature :
                void setColourMap(RDKit::MolDraw2DUtils::ContourParams {lvalue},boost::python::api::object)
        """
    def setContourColour(self, colour: tuple) -> None: 
        """
        setContourColour( self: ContourParams, colour: tuple) -> None

            C++ signature :
                void setContourColour(RDKit::MolDraw2DUtils::ContourParams {lvalue},boost::python::tuple)
        """
    @property
    def contourWidth(self) -> None:
        """
        line width of the contours

        :type: None
        """
    @property
    def dashNegative(self) -> None:
        """
        use a dashed line for negative contours

        :type: None
        """
    @property
    def extraGridPadding(self) -> None:
        """
        extra space (in molecule coords) around the grid

        :type: None
        """
    @property
    def fillGrid(self) -> None:
        """
        colors the grid in addition to drawing contours

        :type: None
        """
    @property
    def gridResolution(self) -> None:
        """
        set the resolution of the grid

        :type: None
        """
    @property
    def setScale(self) -> None:
        """
        set the scale of the drawing object (useful if you draw the grid/contours first)

        :type: None
        """
    __instance_size__ = 112
    pass
class IntStringMap(Boost.Python.instance):
    @staticmethod
    def __contains__( arg1: IntStringMap, arg2: object) -> bool: 
        """
        __contains__( arg1: IntStringMap, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: IntStringMap, arg2: object) -> None: 
        """
        __delitem__( arg1: IntStringMap, arg2: object) -> None

            C++ signature :
                void __delitem__(std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >&>,_object*)
        """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::_Rb_tree_iterator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > __iter__(boost::python::back_reference<std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > >&>)
        """
    @staticmethod
    def __len__( arg1: IntStringMap) -> int: 
        """
        __len__( arg1: IntStringMap) -> int

            C++ signature :
                unsigned long __len__(std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: IntStringMap, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: IntStringMap, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::map<int, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::less<int>, std::allocator<std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > > > {lvalue},_object*,_object*)
        """
    __instance_size__ = 72
    pass
class MolDraw2D(Boost.Python.instance):
    """
    Drawer abstract base class
    """
    def ClearDrawing(self) -> None: 
        """
        ClearDrawing( self: MolDraw2D) -> None
            clears the drawing by filling it with the background color

            C++ signature :
                void ClearDrawing(RDKit::MolDraw2D {lvalue})
        """
    def DrawArc(self, center: Point2D, radius: float, angle1: float, angle2: float, rawCoords: bool = False) -> None: 
        """
        DrawArc( self: MolDraw2D, center: Point2D, radius: float, angle1: float, angle2: float, rawCoords: bool = False) -> None
            draws an arc with the current drawing style. The coordinates are in the molecule frame, the angles are in degrees, angle2 should be > angle1.

            C++ signature :
                void DrawArc(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,double,double,double [,bool=False])
        """
    def DrawArrow(self, cds1: Point2D, cds2: Point2D, asPolygon: bool = False, frac: float = 0.05, angle: float = 0.5235987755982988, color: AtomPairsParameters = None, rawCoords: bool = False) -> None: 
        """
        DrawArrow( self: MolDraw2D, cds1: Point2D, cds2: Point2D, asPolygon: bool = False, frac: float = 0.05, angle: float = 0.5235987755982988, color: AtomPairsParameters = None, rawCoords: bool = False) -> None
            draws an arrow with the current drawing style. The coordinates are in the molecule frame. If asPolygon is true the head of the arrow will be drawn as a triangle, otherwise two lines are used.

            C++ signature :
                void DrawArrow(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False [,double=0.05 [,double=0.5235987755982988 [,boost::python::api::object=None [,bool=False]]]]])
        """
    def DrawAttachmentLine(self, cds1: Point2D, cds2: Point2D, color: tuple, len: float = 1.0, nSegments: int = 16, rawCoords: bool = False) -> None: 
        """
        DrawAttachmentLine( self: MolDraw2D, cds1: Point2D, cds2: Point2D, color: tuple, len: float = 1.0, nSegments: int = 16, rawCoords: bool = False) -> None
            draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond)

            C++ signature :
                void DrawAttachmentLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,boost::python::tuple {lvalue} [,double=1.0 [,unsigned int=16 [,bool=False]]])
        """
    def DrawEllipse(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None: 
        """
        DrawEllipse( self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None
            draws a triangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame

            C++ signature :
                void DrawEllipse(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
    def DrawLine(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None: 
        """
        DrawLine( self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None
            draws a line with the current drawing style. The coordinates are in the molecule frame

            C++ signature :
                void DrawLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
    @typing.overload
    def DrawMolecule(self, mol: Mol, highlightAtoms: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, legend: str = '') -> None: 
        """
        DrawMolecule( self: MolDraw2D, mol: Mol, highlightAtoms: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, legend: str = '') -> None
            renders a molecule
            

            C++ signature :
                void DrawMolecule(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']]]]])

            C++ signature :
                void DrawMolecule(RDKit::MolDraw2D {lvalue},RDKit::ROMol,boost::python::api::object,boost::python::api::object [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']]]]])
        """
    @typing.overload
    def DrawMolecule(self, mol: Mol, highlightAtoms: AtomPairsParameters, highlightBonds: AtomPairsParameters, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, legend: str = '') -> None: ...
    def DrawMoleculeWithHighlights(self, mol: Mol, legend: str, highlight_atom_map: AtomPairsParameters, highlight_bond_map: AtomPairsParameters, highlight_radii: AtomPairsParameters, highlight_linewidth_multipliers: AtomPairsParameters, confId: int = -1) -> None: 
        """
        DrawMoleculeWithHighlights( self: MolDraw2D, mol: Mol, legend: str, highlight_atom_map: AtomPairsParameters, highlight_bond_map: AtomPairsParameters, highlight_radii: AtomPairsParameters, highlight_linewidth_multipliers: AtomPairsParameters, confId: int = -1) -> None
            renders a molecule with multiple highlight colours
            

            C++ signature :
                void DrawMoleculeWithHighlights(RDKit::MolDraw2D {lvalue},RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,boost::python::api::object,boost::python::api::object,boost::python::api::object,boost::python::api::object [,int=-1])
        """
    def DrawMolecules(self, mols: AtomPairsParameters, highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confIds: AtomPairsParameters = None, legends: AtomPairsParameters = None) -> None: 
        """
        DrawMolecules( self: MolDraw2D, mols: AtomPairsParameters, highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confIds: AtomPairsParameters = None, legends: AtomPairsParameters = None) -> None
            renders multiple molecules
            

            C++ signature :
                void DrawMolecules(RDKit::MolDraw2D {lvalue},boost::python::api::object [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]])
        """
    def DrawPolygon(self, cds: AtomPairsParameters, rawCoords: bool = False) -> None: 
        """
        DrawPolygon( self: MolDraw2D, cds: AtomPairsParameters, rawCoords: bool = False) -> None
            draws a polygon with the current drawing style. The coordinates are in the molecule frame

            C++ signature :
                void DrawPolygon(RDKit::MolDraw2D {lvalue},boost::python::api::object [,bool=False])
        """
    def DrawReaction(self, rxn: ChemicalReaction, highlightByReactant: bool = False, highlightColorsReactants: AtomPairsParameters = None, confIds: AtomPairsParameters = None) -> None: 
        """
        DrawReaction( self: MolDraw2D, rxn: ChemicalReaction, highlightByReactant: bool = False, highlightColorsReactants: AtomPairsParameters = None, confIds: AtomPairsParameters = None) -> None
            renders a reaction
            

            C++ signature :
                void DrawReaction(RDKit::MolDraw2D {lvalue},RDKit::ChemicalReaction [,bool=False [,boost::python::api::object=None [,boost::python::api::object=None]]])
        """
    def DrawRect(self, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None: 
        """
        DrawRect( self: MolDraw2D, cds1: Point2D, cds2: Point2D, rawCoords: bool = False) -> None
            draws a rectangle with the current drawing style in the rectangle defined by the two points. The coordinates are in the molecule frame

            C++ signature :
                void DrawRect(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
    @typing.overload
    def DrawString(self, string: str, pos: Point2D, rawCoords: bool = False) -> None: 
        """
        DrawString( self: MolDraw2D, string: str, pos: Point2D, rawCoords: bool = False) -> None
            add text to the canvas

            C++ signature :
                void DrawString(RDKit::MolDraw2D {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point2D [,bool=False])

            C++ signature :
                void DrawString(RDKit::MolDraw2D {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point2D,int [,bool=False])
        """
    @typing.overload
    def DrawString(self, string: str, pos: Point2D, align: int, rawCoords: bool = False) -> None: ...
    def DrawTriangle(self, cds1: Point2D, cds2: Point2D, cds3: Point2D, rawCoords: bool = False) -> None: 
        """
        DrawTriangle( self: MolDraw2D, cds1: Point2D, cds2: Point2D, cds3: Point2D, rawCoords: bool = False) -> None
            draws a triangle with the current drawing style. The coordinates are in the molecule frame

            C++ signature :
                void DrawTriangle(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,RDGeom::Point2D [,bool=False])
        """
    def DrawWavyLine(self, cds1: Point2D, cds2: Point2D, color1: tuple, color2: tuple, nSegments: int = 16, vertOffset: float = 0.05, rawCoords: bool = False) -> None: 
        """
        DrawWavyLine( self: MolDraw2D, cds1: Point2D, cds2: Point2D, color1: tuple, color2: tuple, nSegments: int = 16, vertOffset: float = 0.05, rawCoords: bool = False) -> None
            draw a line indicating the presence of an attachment point (normally a squiggle line perpendicular to a bond)

            C++ signature :
                void DrawWavyLine(RDKit::MolDraw2D {lvalue},RDGeom::Point2D,RDGeom::Point2D,boost::python::tuple {lvalue},boost::python::tuple {lvalue} [,unsigned int=16 [,double=0.05 [,bool=False]]])
        """
    @staticmethod
    def FillPolys( arg1: MolDraw2D) -> bool: 
        """
        FillPolys( arg1: MolDraw2D) -> bool
            returns whether or not polygons are being filled

            C++ signature :
                bool FillPolys(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def FlexiMode( arg1: MolDraw2D) -> bool: 
        """
        FlexiMode( arg1: MolDraw2D) -> bool
            returns whether or not FlexiMode is being used

            C++ signature :
                bool FlexiMode(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def FontSize( arg1: MolDraw2D) -> float: 
        """
        FontSize( arg1: MolDraw2D) -> float
            get the default font size. The units are, roughly, pixels.

            C++ signature :
                double FontSize(RDKit::MolDraw2D {lvalue})
        """
    @typing.overload
    def GetDrawCoords(self, point: Point2D) -> Point2D: 
        """
        GetDrawCoords( self: MolDraw2D, point: Point2D) -> Point2D
            get the coordinates in drawing space for a particular point in molecule space

            C++ signature :
                RDGeom::Point2D GetDrawCoords(RDKit::MolDraw2D {lvalue},RDGeom::Point2D)

            C++ signature :
                RDGeom::Point2D GetDrawCoords(RDKit::MolDraw2D {lvalue},int)
        """
    @typing.overload
    def GetDrawCoords(self, atomIndex: int) -> Point2D: ...
    def GetMolSize(self, mol: Mol, highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, legend: str = '') -> tuple: 
        """
        GetMolSize( self: MolDraw2D, mol: Mol, highlightAtoms: AtomPairsParameters = None, highlightBonds: AtomPairsParameters = None, highlightAtomColors: AtomPairsParameters = None, highlightBondColors: AtomPairsParameters = None, highlightAtomRadii: AtomPairsParameters = None, confId: int = -1, legend: str = '') -> tuple
            returns the width and height required to draw a molecule at the current size

            C++ signature :
                boost::python::tuple GetMolSize(RDKit::MolDraw2D {lvalue},RDKit::ROMol [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::api::object=None [,int=-1 [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='']]]]]]])
        """
    @staticmethod
    def Height( arg1: MolDraw2D) -> int: 
        """
        Height( arg1: MolDraw2D) -> int
            get the height of the drawing canvas

            C++ signature :
                int Height(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def LineWidth( arg1: MolDraw2D) -> float: 
        """
        LineWidth( arg1: MolDraw2D) -> float
            returns the line width being used

            C++ signature :
                double LineWidth(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def Offset( arg1: MolDraw2D) -> Point2D: 
        """
        Offset( arg1: MolDraw2D) -> Point2D
            returns the offset (in drawing coordinates) for the drawing

            C++ signature :
                RDGeom::Point2D Offset(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def SetColour( arg1: MolDraw2D, arg2: tuple) -> None: 
        """
        SetColour( arg1: MolDraw2D, arg2: tuple) -> None
            set the color being used fr drawing and filling

            C++ signature :
                void SetColour(RDKit::MolDraw2D {lvalue},boost::python::tuple)
        """
    @staticmethod
    def SetDrawOptions( arg1: MolDraw2D, arg2: MolDrawOptions) -> None: 
        """
        SetDrawOptions( arg1: MolDraw2D, arg2: MolDrawOptions) -> None
            Copies the drawing options passed in over our drawing options

            C++ signature :
                void SetDrawOptions(RDKit::MolDraw2D {lvalue},RDKit::MolDrawOptions)
        """
    @staticmethod
    def SetFillPolys( arg1: MolDraw2D, arg2: bool) -> None: 
        """
        SetFillPolys( arg1: MolDraw2D, arg2: bool) -> None
            sets whether or not polygons are filled

            C++ signature :
                void SetFillPolys(RDKit::MolDraw2D {lvalue},bool)
        """
    @staticmethod
    def SetFlexiMode( arg1: MolDraw2D, arg2: bool) -> None: 
        """
        SetFlexiMode( arg1: MolDraw2D, arg2: bool) -> None
            when FlexiMode is set, molecules will always been drawn with the default values for bond length, font size, etc.

            C++ signature :
                void SetFlexiMode(RDKit::MolDraw2D {lvalue},bool)
        """
    @staticmethod
    def SetFontSize( arg1: MolDraw2D, arg2: float) -> None: 
        """
        SetFontSize( arg1: MolDraw2D, arg2: float) -> None
            change the default font size. The units are, roughly, pixels.

            C++ signature :
                void SetFontSize(RDKit::MolDraw2D {lvalue},double)
        """
    @staticmethod
    def SetLineWidth( arg1: MolDraw2D, arg2: float) -> None: 
        """
        SetLineWidth( arg1: MolDraw2D, arg2: float) -> None
            set the line width being used

            C++ signature :
                void SetLineWidth(RDKit::MolDraw2D {lvalue},double)
        """
    @staticmethod
    def SetOffset( arg1: MolDraw2D, arg2: int, arg3: int) -> None: 
        """
        SetOffset( arg1: MolDraw2D, arg2: int, arg3: int) -> None
            set the offset (in drawing coordinates) for the drawing

            C++ signature :
                void SetOffset(RDKit::MolDraw2D {lvalue},int,int)
        """
    def SetScale(self, width: int, height: int, minv: Point2D, maxv: Point2D, mol: AtomPairsParameters = None) -> None: 
        """
        SetScale( self: MolDraw2D, width: int, height: int, minv: Point2D, maxv: Point2D, mol: AtomPairsParameters = None) -> None
            uses the values provided to set the drawing scaling

            C++ signature :
                void SetScale(RDKit::MolDraw2D {lvalue},int,int,RDGeom::Point2D,RDGeom::Point2D [,boost::python::api::object=None])
        """
    @staticmethod
    def Width( arg1: MolDraw2D) -> int: 
        """
        Width( arg1: MolDraw2D) -> int
            get the width of the drawing canvas

            C++ signature :
                int Width(RDKit::MolDraw2D {lvalue})
        """
    @staticmethod
    def drawOptions( arg1: MolDraw2D) -> MolDrawOptions: 
        """
        drawOptions( arg1: MolDraw2D) -> MolDrawOptions
            Returns a modifiable version of the current drawing options

            C++ signature :
                RDKit::MolDrawOptions {lvalue} drawOptions(RDKit::MolDraw2D {lvalue})
        """
    pass
class MolDraw2DCairo(MolDraw2D, Boost.Python.instance):
    """
    Cairo molecule drawer
    """
    @staticmethod
    def FinishDrawing( arg1: MolDraw2DCairo) -> None: 
        """
        FinishDrawing( arg1: MolDraw2DCairo) -> None
            add the last bits to finish the drawing

            C++ signature :
                void FinishDrawing(RDKit::MolDraw2DCairo {lvalue})
        """
    @staticmethod
    def GetDrawingText( arg1: MolDraw2DCairo) -> object: 
        """
        GetDrawingText( arg1: MolDraw2DCairo) -> object
            return the PNG data as a string

            C++ signature :
                boost::python::api::object GetDrawingText(RDKit::MolDraw2DCairo)
        """
    @staticmethod
    def WriteDrawingText( arg1: MolDraw2DCairo, arg2: str) -> None: 
        """
        WriteDrawingText( arg1: MolDraw2DCairo, arg2: str) -> None
            write the PNG data to the named file

            C++ signature :
                void WriteDrawingText(RDKit::MolDraw2DCairo {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
    @staticmethod
    def __init__( arg1: object, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None: 
        """
        __init__( arg1: object, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None

            C++ signature :
                void __init__(_object*,int,int [,int=-1 [,int=-1 [,bool=False]]])
        """
    __instance_size__ = 904
    pass
class MolDraw2DSVG(MolDraw2D, Boost.Python.instance):
    """
    SVG molecule drawer
    """
    @staticmethod
    def AddMoleculeMetadata( arg1: MolDraw2DSVG, mol: Mol, confId: int = -1) -> None: 
        """
        AddMoleculeMetadata( arg1: MolDraw2DSVG, mol: Mol, confId: int = -1) -> None
            add RDKit-specific information to the bottom of the drawing

            C++ signature :
                void AddMoleculeMetadata(RDKit::MolDraw2DSVG {lvalue},RDKit::ROMol [,int=-1])
        """
    @staticmethod
    def FinishDrawing( arg1: MolDraw2DSVG) -> None: 
        """
        FinishDrawing( arg1: MolDraw2DSVG) -> None
            add the last bits of SVG to finish the drawing

            C++ signature :
                void FinishDrawing(RDKit::MolDraw2DSVG {lvalue})
        """
    @staticmethod
    def GetDrawingText( arg1: MolDraw2DSVG) -> str: 
        """
        GetDrawingText( arg1: MolDraw2DSVG) -> str
            return the SVG

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetDrawingText(RDKit::MolDraw2DSVG {lvalue})
        """
    @staticmethod
    def TagAtoms( arg1: MolDraw2DSVG, mol: Mol, radius: float = 0.2, events: AtomPairsParameters = None) -> None: 
        """
        TagAtoms( arg1: MolDraw2DSVG, mol: Mol, radius: float = 0.2, events: AtomPairsParameters = None) -> None
            allow atom selection in the SVG

            C++ signature :
                void TagAtoms(RDKit::MolDraw2DSVG {lvalue},RDKit::ROMol [,double=0.2 [,boost::python::api::object=None]])
        """
    @staticmethod
    def __init__( arg1: object, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None: 
        """
        __init__( arg1: object, width: int, height: int, panelWidth: int = -1, panelHeight: int = -1, noFreetype: bool = False) -> None

            C++ signature :
                void __init__(_object*,int,int [,int=-1 [,int=-1 [,bool=False]]])
        """
    __instance_size__ = 1288
    pass
class MolDrawOptions(Boost.Python.instance):
    """
    Drawing options
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def getAnnotationColour( arg1: MolDrawOptions) -> object: 
        """
        getAnnotationColour( arg1: MolDrawOptions) -> object
            method returning the annotation colour

            C++ signature :
                boost::python::api::object getAnnotationColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getBackgroundColour( arg1: MolDrawOptions) -> object: 
        """
        getBackgroundColour( arg1: MolDrawOptions) -> object
            method returning the background colour

            C++ signature :
                boost::python::api::object getBackgroundColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getHighlightColour( arg1: MolDrawOptions) -> object: 
        """
        getHighlightColour( arg1: MolDrawOptions) -> object
            method returning the highlight colour

            C++ signature :
                boost::python::api::object getHighlightColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getLegendColour( arg1: MolDrawOptions) -> object: 
        """
        getLegendColour( arg1: MolDrawOptions) -> object
            method returning the legend colour

            C++ signature :
                boost::python::api::object getLegendColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getQueryColour( arg1: MolDrawOptions) -> object: 
        """
        getQueryColour( arg1: MolDrawOptions) -> object
            method returning the query colour

            C++ signature :
                boost::python::api::object getQueryColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getSymbolColour( arg1: MolDrawOptions) -> object: 
        """
        getSymbolColour( arg1: MolDrawOptions) -> object
            method returning the symbol colour

            C++ signature :
                boost::python::api::object getSymbolColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def getVariableAttachmentColour( arg1: MolDrawOptions) -> object: 
        """
        getVariableAttachmentColour( arg1: MolDrawOptions) -> object
            method for getting the colour of variable attachment points

            C++ signature :
                boost::python::api::object getVariableAttachmentColour(RDKit::MolDrawOptions)
        """
    @staticmethod
    def setAnnotationColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setAnnotationColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the annotation colour

            C++ signature :
                void setAnnotationColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setAtomPalette( arg1: MolDrawOptions, arg2: AtomPairsParameters) -> None: 
        """
        setAtomPalette( arg1: MolDrawOptions, arg2: AtomPairsParameters) -> None
            sets the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

            C++ signature :
                void setAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)
        """
    @staticmethod
    def setBackgroundColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setBackgroundColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the background colour

            C++ signature :
                void setBackgroundColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setHighlightColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setHighlightColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the highlight colour

            C++ signature :
                void setHighlightColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setLegendColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setLegendColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the legend colour

            C++ signature :
                void setLegendColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setQueryColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setQueryColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the query colour

            C++ signature :
                void setQueryColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setSymbolColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setSymbolColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the symbol colour

            C++ signature :
                void setSymbolColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def setVariableAttachmentColour( arg1: MolDrawOptions, arg2: tuple) -> None: 
        """
        setVariableAttachmentColour( arg1: MolDrawOptions, arg2: tuple) -> None
            method for setting the colour of variable attachment points

            C++ signature :
                void setVariableAttachmentColour(RDKit::MolDrawOptions {lvalue},boost::python::tuple)
        """
    @staticmethod
    def updateAtomPalette( arg1: MolDrawOptions, arg2: AtomPairsParameters) -> None: 
        """
        updateAtomPalette( arg1: MolDrawOptions, arg2: AtomPairsParameters) -> None
            updates the palette for atoms and bonds from a dictionary mapping ints to 3-tuples

            C++ signature :
                void updateAtomPalette(RDKit::MolDrawOptions {lvalue},boost::python::api::object)
        """
    @staticmethod
    def useAvalonAtomPalette( arg1: MolDrawOptions) -> None: 
        """
        useAvalonAtomPalette( arg1: MolDrawOptions) -> None
            use the Avalon renderer palette for atoms and bonds

            C++ signature :
                void useAvalonAtomPalette(RDKit::MolDrawOptions {lvalue})
        """
    @staticmethod
    def useBWAtomPalette( arg1: MolDrawOptions) -> None: 
        """
        useBWAtomPalette( arg1: MolDrawOptions) -> None
            use a black and white palette for atoms and bonds

            C++ signature :
                void useBWAtomPalette(RDKit::MolDrawOptions {lvalue})
        """
    @staticmethod
    def useCDKAtomPalette( arg1: MolDrawOptions) -> None: 
        """
        useCDKAtomPalette( arg1: MolDrawOptions) -> None
            use the CDK palette for atoms and bonds

            C++ signature :
                void useCDKAtomPalette(RDKit::MolDrawOptions {lvalue})
        """
    @staticmethod
    def useDefaultAtomPalette( arg1: MolDrawOptions) -> None: 
        """
        useDefaultAtomPalette( arg1: MolDrawOptions) -> None
            use the default colour palette for atoms and bonds

            C++ signature :
                void useDefaultAtomPalette(RDKit::MolDrawOptions {lvalue})
        """
    @property
    def addAtomIndices(self) -> None:
        """
        adds atom indices to drawings. Default False.

        :type: None
        """
    @property
    def addBondIndices(self) -> None:
        """
        adds bond indices to drawings. Default False.

        :type: None
        """
    @property
    def addStereoAnnotation(self) -> None:
        """
        adds R/S and E/Z to drawings. Default False.

        :type: None
        """
    @property
    def additionalAtomLabelPadding(self) -> None:
        """
        additional padding to leave around atom labels. Expressed as a fraction of the font size.

        :type: None
        """
    @property
    def annotationFontScale(self) -> None:
        """
        Scale of font for atom and bond annotation relative to atomlabel font.  Default=0.75.

        :type: None
        """
    @property
    def atomHighlightsAreCircles(self) -> None:
        """
        forces atom highlights always to be circles.Default (false) is to put ellipses roundlonger labels.

        :type: None
        """
    @property
    def atomLabelDeuteriumTritium(self) -> None:
        """
        labels deuterium as D and tritium as T

        :type: None
        """
    @property
    def atomLabels(self) -> None:
        """
        maps indices to atom labels

        :type: None
        """
    @property
    def atomRegions(self) -> None:
        """
        regions to outline

        :type: None
        """
    @property
    def baseFontSize(self) -> None:
        """
        relative size of font.  Defaults to 0.6.  -1 means use default.

        :type: None
        """
    @property
    def bondLineWidth(self) -> None:
        """
        if positive, this overrides the default line width for bonds

        :type: None
        """
    @property
    def centreMoleculesBeforeDrawing(self) -> None:
        """
        Moves the centre of the drawn molecule to (0,0).Default False.

        :type: None
        """
    @property
    def circleAtoms(self) -> None:
        """
        :type: None
        """
    @property
    def clearBackground(self) -> None:
        """
        clear the background before drawing a molecule

        :type: None
        """
    @property
    def comicMode(self) -> None:
        """
        simulate hand-drawn lines for bonds. When combined with a font like Comic-Sans or Comic-Neue, this gives xkcd-like drawings. Default is false.

        :type: None
        """
    @property
    def continuousHighlight(self) -> None:
        """
        :type: None
        """
    @property
    def drawMolsSameScale(self) -> None:
        """
        when drawing multiple molecules with DrawMolecules, forces them to use the same scale.  Default is true.

        :type: None
        """
    @property
    def dummiesAreAttachments(self) -> None:
        """
        :type: None
        """
    @property
    def dummyIsotopeLabels(self) -> None:
        """
        adds isotope labels on dummy atoms. Default True.

        :type: None
        """
    @property
    def explicitMethyl(self) -> None:
        """
        Draw terminal methyls explictly.  Default is false.

        :type: None
        """
    @property
    def fillHighlights(self) -> None:
        """
        :type: None
        """
    @property
    def fixedBondLength(self) -> None:
        """
        If > 0.0, fixes bond length to this number of pixelsunless that would make it too big.  Default -1.0 meansno fix.  If both set, fixedScale takes precedence.

        :type: None
        """
    @property
    def fixedFontSize(self) -> None:
        """
        font size in pixels. default=-1 means not fixed.  If set, always used irrespective of scale, minFontSize and maxFontSize.

        :type: None
        """
    @property
    def fixedScale(self) -> None:
        """
        If > 0.0, fixes scale to that fraction of width ofdraw window.  Default -1.0 means adjust scale to fit.

        :type: None
        """
    @property
    def flagCloseContactsDist(self) -> None:
        """
        :type: None
        """
    @property
    def fontFile(self) -> None:
        """
        Font file for use with FreeType text drawer.  Can also be BuiltinTelexRegular (the default) or BuiltinRobotoRegular.

        :type: None
        """
    @property
    def highlightBondWidthMultiplier(self) -> None:
        """
        What to multiply default bond width by for highlighting bonds. Default-8.

        :type: None
        """
    @property
    def highlightRadius(self) -> None:
        """
        Default radius for highlight circles.

        :type: None
        """
    @property
    def includeAtomTags(self) -> None:
        """
        include atom tags in output

        :type: None
        """
    @property
    def includeChiralFlagLabel(self) -> None:
        """
        add a molecule annotation with "ABS" if the chiral flag is set. Default is false.

        :type: None
        """
    @property
    def includeMetadata(self) -> None:
        """
        When possible, include metadata about molecules and reactions to allow them to be reconstructed. Default is true.

        :type: None
        """
    @property
    def includeRadicals(self) -> None:
        """
        include radicals in the drawing (it can be useful to turn this off for reactions and queries). Default is true.

        :type: None
        """
    @property
    def isotopeLabels(self) -> None:
        """
        adds isotope labels on non-dummy atoms. Default True.

        :type: None
        """
    @property
    def legendFontSize(self) -> None:
        """
        font size in pixels of the legend (if drawn)

        :type: None
        """
    @property
    def legendFraction(self) -> None:
        """
        fraction of the draw panel to be used for the legend if present

        :type: None
        """
    @property
    def maxFontSize(self) -> None:
        """
        maximum font size in pixels. default=40, -1 means no maximum.

        :type: None
        """
    @property
    def minFontSize(self) -> None:
        """
        minimum font size in pixels. default=6, -1 means no minimum.

        :type: None
        """
    @property
    def multiColourHighlightStyle(self) -> None:
        """
        Either 'CircleAndLine' or 'Lasso', to control style ofmulti-coloured highlighting in DrawMoleculeWithHighlights.Default is CircleAndLine.

        :type: None
        """
    @property
    def multipleBondOffset(self) -> None:
        """
        offset for the extra lines in a multiple bond as a fraction of mean bond length

        :type: None
        """
    @property
    def noAtomLabels(self) -> None:
        """
        disables inclusion of atom labels in the rendering

        :type: None
        """
    @property
    def padding(self) -> None:
        """
        fraction of empty space to leave around molecule

        :type: None
        """
    @property
    def prepareMolsBeforeDrawing(self) -> None:
        """
        call prepareMolForDrawing() on each molecule passed to DrawMolecules()

        :type: None
        """
    @property
    def rotate(self) -> None:
        """
        Rotates molecule about centre by this number of degrees,

        :type: None
        """
    @property
    def scaleBondWidth(self) -> None:
        """
        Scales the width of drawn bonds using image scaling.

        :type: None
        """
    @property
    def scaleHighlightBondWidth(self) -> None:
        """
        Scales the width of drawn highlighted bonds using image scaling.

        :type: None
        """
    @property
    def scalingFactor(self) -> None:
        """
        scaling factor for pixels->angstrom when auto scalingbeing used.  Default is 20.

        :type: None
        """
    @property
    def simplifiedStereoGroupLabel(self) -> None:
        """
        if all specified stereocenters are in a single StereoGroup, show a molecule-level annotation instead of the individual labels. Default is false.

        :type: None
        """
    @property
    def singleColourWedgeBonds(self) -> None:
        """
        if true wedged and dashed bonds are drawn using symbolColour rather than inheriting their colour from the atoms. Default is false.

        :type: None
        """
    @property
    def splitBonds(self) -> None:
        """
        :type: None
        """
    @property
    def unspecifiedStereoIsUnknown(self) -> None:
        """
        if true, double bonds with unspecified stereo are drawn crossed, potential stereocenters with unspecified stereo are drawn with a wavy bond. Default is false.

        :type: None
        """
    @property
    def useComplexQueryAtomSymbols(self) -> None:
        """
        replace any atom, any hetero, any halo queries with complex query symbols A, Q, X, M, optionally followed by H if hydrogen is included (except for AH, which stays *). Default is true

        :type: None
        """
    @property
    def useMolBlockWedging(self) -> None:
        """
        If the molecule came from a MolBlock, prefer the wedging information that provides.  If false, use RDKit rules.  Default false

        :type: None
        """
    @property
    def variableAtomRadius(self) -> None:
        """
        radius value to use for atoms involved in variable attachment points.

        :type: None
        """
    @property
    def variableBondWidthMultiplier(self) -> None:
        """
        what to multiply standard bond width by for variable attachment points.

        :type: None
        """
    __instance_size__ = 640
    pass
class MultiColourHighlightStyle(Boost.Python.enum, int):
    CircleAndLine = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine
    Lasso = rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso
    __slots__ = ()
    names = {'CircleAndLine': rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine, 'Lasso': rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso}
    values = {0: rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.CircleAndLine, 1: rdkit.Chem.Draw.rdMolDraw2D.MultiColourHighlightStyle.Lasso}
    pass
class map_indexing_suite_IntStringMap_entry(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __repr__( arg1: map_indexing_suite_IntStringMap_entry) -> object: 
        """
        __repr__( arg1: map_indexing_suite_IntStringMap_entry) -> object

            C++ signature :
                boost::python::api::object __repr__(std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >)
        """
    @staticmethod
    def data( arg1: map_indexing_suite_IntStringMap_entry) -> str: 
        """
        data( arg1: map_indexing_suite_IntStringMap_entry) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > data(std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > {lvalue})
        """
    @staticmethod
    def key( arg1: map_indexing_suite_IntStringMap_entry) -> int: 
        """
        key( arg1: map_indexing_suite_IntStringMap_entry) -> int

            C++ signature :
                int key(std::pair<int const, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > {lvalue})
        """
    __instance_size__ = 64
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
            void ContourAndDrawGaussians(RDKit::MolDraw2D {lvalue},boost::python::api::object,boost::python::api::object,boost::python::api::object [,unsigned int=10 [,boost::python::api::object=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7fd05bef1640> [,boost::python::api::object=None]]]])
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
            void ContourAndDrawGrid(RDKit::MolDraw2D {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue},boost::python::api::object {lvalue} [,unsigned int=10 [,boost::python::api::object {lvalue}=None [,RDKit::MolDraw2DUtils::ContourParams=<rdkit.Chem.Draw.rdMolDraw2D.ContourParams object at 0x7fd05bef1700> [,boost::python::api::object=None]]]])
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
