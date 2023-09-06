"""Module containing tools for normalizing molecules defined by SMARTS patterns"""
from __future__ import annotations
import rdkit.Chem.MolStandardize.rdMolStandardize
import typing
import Boost.Python

__all__ = [
    "AllowedAtomsValidation",
    "CHARGE_CORRECTIONS",
    "CanonicalTautomer",
    "ChargeCorrection",
    "ChargeParent",
    "Cleanup",
    "CleanupInPlace",
    "CleanupParameters",
    "DisallowedAtomsValidation",
    "DisconnectOrganometallics",
    "DisconnectOrganometallicsInPlace",
    "FragmentParent",
    "FragmentRemover",
    "FragmentRemoverFromData",
    "FragmentValidation",
    "GetV1TautomerEnumerator",
    "IsotopeParent",
    "IsotopeValidation",
    "LargestFragmentChooser",
    "MetalDisconnector",
    "MetalDisconnectorOptions",
    "MolVSValidation",
    "MolVSValidations",
    "NeutralValidation",
    "NoAtomValidation",
    "Normalize",
    "NormalizeInPlace",
    "Normalizer",
    "NormalizerFromData",
    "NormalizerFromParams",
    "RDKitValidation",
    "Reionize",
    "ReionizeInPlace",
    "Reionizer",
    "ReionizerFromData",
    "RemoveFragments",
    "RemoveFragmentsInPlace",
    "SmilesTautomerMap",
    "StandardizeSmiles",
    "StereoParent",
    "SuperParent",
    "Tautomer",
    "TautomerEnumerator",
    "TautomerEnumeratorCallback",
    "TautomerEnumeratorResult",
    "TautomerEnumeratorStatus",
    "TautomerParent",
    "Uncharger",
    "UpdateParamsFromJSON",
    "ValidateSmiles",
    "map_indexing_suite_SmilesTautomerMap_entry"
]


class AllowedAtomsValidation(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: AtomPairsParameters, arg2: AtomPairsParameters) -> object: 
        """
        __init__( arg1: AtomPairsParameters, arg2: AtomPairsParameters) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object)
        """
    def validate(self, mol: Mol, reportAllFailures: bool = False) -> list: 
        """
        validate( self: AllowedAtomsValidation, mol: Mol, reportAllFailures: bool = False) -> list

            C++ signature :
                boost::python::list validate(RDKit::MolStandardize::AllowedAtomsValidation {lvalue},RDKit::ROMol [,bool=False])
        """
    pass
class ChargeCorrection(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object, name: str, smarts: str, charge: int) -> None: 
        """
        __init__( arg1: object, name: str, smarts: str, charge: int) -> None

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int)
        """
    @property
    def Charge(self) -> None:
        """
        :type: None
        """
    @property
    def Name(self) -> None:
        """
        :type: None
        """
    @property
    def Smarts(self) -> None:
        """
        :type: None
        """
    __instance_size__ = 96
    pass
class CleanupParameters(Boost.Python.instance):
    """
    Parameters controlling molecular standardization
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def acidbaseFile(self) -> None:
        """
        file containing the acid and base definitions

        :type: None
        """
    @property
    def doCanonical(self) -> None:
        """
        apply atom-order dependent normalizations (like uncharging) in a canonical order

        :type: None
        """
    @property
    def fragmentFile(self) -> None:
        """
        file containing the acid and base definitions

        :type: None
        """
    @property
    def largestFragmentChooserCountHeavyAtomsOnly(self) -> None:
        """
        whether LargestFragmentChooser should only count heavy atoms (defaults to False)

        :type: None
        """
    @property
    def largestFragmentChooserUseAtomCount(self) -> None:
        """
        Whether LargestFragmentChooser should use atom count as main criterion before MW (defaults to True)

        :type: None
        """
    @property
    def maxRestarts(self) -> None:
        """
        maximum number of restarts

        :type: None
        """
    @property
    def maxTautomers(self) -> None:
        """
        maximum number of tautomers to generate (defaults to 1000)

        :type: None
        """
    @property
    def maxTransforms(self) -> None:
        """
        maximum number of transforms to apply during tautomer enumeration (defaults to 1000)

        :type: None
        """
    @property
    def normalizationsFile(self) -> None:
        """
        file containing the normalization transformations

        :type: None
        """
    @property
    def preferOrganic(self) -> None:
        """
        prefer organic fragments to inorganic ones when deciding what to keep

        :type: None
        """
    @property
    def tautomerReassignStereo(self) -> None:
        """
        call AssignStereochemistry on all generated tautomers (defaults to True)

        :type: None
        """
    @property
    def tautomerRemoveBondStereo(self) -> None:
        """
        remove stereochemistry from double bonds involved in tautomerism (defaults to True)

        :type: None
        """
    @property
    def tautomerRemoveIsotopicHs(self) -> None:
        """
        remove isotopic Hs from centers involved in tautomerism (defaults to True)

        :type: None
        """
    @property
    def tautomerRemoveSp3Stereo(self) -> None:
        """
        remove stereochemistry from sp3 centers involved in tautomerism (defaults to True)

        :type: None
        """
    @property
    def tautomerTransformsFile(self) -> None:
        """
        file containing the tautomer transformations

        :type: None
        """
    __instance_size__ = 312
    pass
class DisallowedAtomsValidation(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: AtomPairsParameters, arg2: AtomPairsParameters) -> object: 
        """
        __init__( arg1: AtomPairsParameters, arg2: AtomPairsParameters) -> object

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object)
        """
    def validate(self, mol: Mol, reportAllFailures: bool = False) -> list: 
        """
        validate( self: DisallowedAtomsValidation, mol: Mol, reportAllFailures: bool = False) -> list

            C++ signature :
                boost::python::list validate(RDKit::MolStandardize::DisallowedAtomsValidation {lvalue},RDKit::ROMol [,bool=False])
        """
    pass
class FragmentRemover(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object* [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >='' [,bool=True [,bool=False]]])
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, fragmentFilename: str = '', leave_last: bool = True, skip_if_all_match: bool = False) -> None: ...
    def remove(self, mol: Mol) -> Mol: 
        """
        remove( self: FragmentRemover, mol: Mol) -> Mol

            C++ signature :
                RDKit::ROMol* remove(RDKit::MolStandardize::FragmentRemover {lvalue},RDKit::ROMol)
        """
    def removeInPlace(self, mol: Mol) -> None: 
        """
        removeInPlace( self: FragmentRemover, mol: Mol) -> None
            modifies the molecule in place

            C++ signature :
                void removeInPlace(RDKit::MolStandardize::FragmentRemover {lvalue},RDKit::ROMol {lvalue})
        """
    __instance_size__ = 40
    pass
class MolVSValidations(Boost.Python.instance):
    def run(self, mol: Mol, reportAllFailures: bool, errors: object) -> None: 
        """
        run( self: MolVSValidations, mol: Mol, reportAllFailures: bool, errors: object) -> None

            C++ signature :
                void run(RDKit::MolStandardize::MolVSValidations {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
    pass
class IsotopeValidation(MolVSValidations, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def run(self, mol: Mol, reportAllFailures: bool, errors: object) -> None: 
        """
        run( self: IsotopeValidation, mol: Mol, reportAllFailures: bool, errors: object) -> None

            C++ signature :
                void run(RDKit::MolStandardize::IsotopeValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
    __instance_size__ = 32
    pass
class LargestFragmentChooser(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object, preferOrganic: bool = False) -> None: 
        """
        __init__( arg1: object, preferOrganic: bool = False) -> None

            C++ signature :
                void __init__(_object* [,bool=False])

            C++ signature :
                void __init__(_object*,RDKit::MolStandardize::CleanupParameters)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, params: CleanupParameters) -> None: ...
    def choose(self, mol: Mol) -> Mol: 
        """
        choose( self: LargestFragmentChooser, mol: Mol) -> Mol

            C++ signature :
                RDKit::ROMol* choose(RDKit::MolStandardize::LargestFragmentChooser {lvalue},RDKit::ROMol)
        """
    __instance_size__ = 32
    pass
class MetalDisconnector(Boost.Python.instance):
    """
    a class to disconnect metals that are defined as covalently bonded to non-metals
    """
    def Disconnect(self, mol: Mol) -> Mol: 
        """
        Disconnect( self: MetalDisconnector, mol: Mol) -> Mol
            performs the disconnection

            C++ signature :
                RDKit::ROMol* Disconnect((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
    def DisconnectInPlace(self, mol: Mol) -> None: 
        """
        DisconnectInPlace( self: MetalDisconnector, mol: Mol) -> None
            performs the disconnection, modifies the input molecule

            C++ signature :
                void DisconnectInPlace((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol {lvalue})
        """
    def SetMetalNof(self, mol: Mol) -> None: 
        """
        SetMetalNof( self: MetalDisconnector, mol: Mol) -> None
            Set the query molecule defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine.

            C++ signature :
                void SetMetalNof((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
    def SetMetalNon(self, mol: Mol) -> None: 
        """
        SetMetalNon( self: MetalDisconnector, mol: Mol) -> None
            Set the query molecule defining the metals to disconnect from other inorganic elements.

            C++ signature :
                void SetMetalNon((anonymous namespace)::MetalDisconnectorWrap {lvalue},RDKit::ROMol)
        """
    @staticmethod
    def __init__( arg1: object, options: AtomPairsParameters = None) -> None: 
        """
        __init__( arg1: object, options: AtomPairsParameters = None) -> None

            C++ signature :
                void __init__(_object* [,boost::python::api::object=None])
        """
    @property
    def MetalNof(self) -> None:
        """
        SMARTS defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine

        :type: None
        """
    @property
    def MetalNon(self) -> None:
        """
        SMARTS defining the metals to disconnect other inorganic elements

        :type: None
        """
    __instance_size__ = 32
    pass
class MetalDisconnectorOptions(Boost.Python.instance):
    """
    Metal Disconnector Options
    """
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @property
    def adjustCharges(self) -> None:
        """
        Whether to adjust charges on ligand atoms.  Default true.

        :type: None
        """
    @property
    def removeHapticDummies(self) -> None:
        """
        Whether to remove the dummy atoms representing haptic bonds.  Such dummies are bonded to the metal with a bond that has the MolFileBondEndPts prop set.  Default false.

        :type: None
        """
    @property
    def splitAromaticC(self) -> None:
        """
        Whether to split metal-aromatic C bonds.  Default false.

        :type: None
        """
    @property
    def splitGrignards(self) -> None:
        """
        Whether to split Grignard-type complexes. Default false.

        :type: None
        """
    __instance_size__ = 32
    pass
class MolVSValidation(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void* __init__(boost::python::api::object,boost::python::api::object)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, arg2: AtomPairsParameters) -> object: ...
    def validate(self, mol: Mol, reportAllFailures: bool = False) -> list: 
        """
        validate( self: MolVSValidation, mol: Mol, reportAllFailures: bool = False) -> list

            C++ signature :
                boost::python::list validate(RDKit::MolStandardize::MolVSValidation {lvalue},RDKit::ROMol [,bool=False])
        """
    __instance_size__ = 56
    pass
class FragmentValidation(MolVSValidations, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def run(self, mol: Mol, reportAllFailures: bool, errors: object) -> None: 
        """
        run( self: FragmentValidation, mol: Mol, reportAllFailures: bool, errors: object) -> None

            C++ signature :
                void run(RDKit::MolStandardize::FragmentValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
    __instance_size__ = 32
    pass
class NeutralValidation(MolVSValidations, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def run(self, mol: Mol, reportAllFailures: bool, errors: object) -> None: 
        """
        run( self: NeutralValidation, mol: Mol, reportAllFailures: bool, errors: object) -> None

            C++ signature :
                void run(RDKit::MolStandardize::NeutralValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
    __instance_size__ = 32
    pass
class NoAtomValidation(MolVSValidations, Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def run(self, mol: Mol, reportAllFailures: bool, errors: object) -> None: 
        """
        run( self: NoAtomValidation, mol: Mol, reportAllFailures: bool, errors: object) -> None

            C++ signature :
                void run(RDKit::MolStandardize::NoAtomValidation {lvalue},RDKit::ROMol,bool,std::vector<RDKit::MolStandardize::ValidationErrorInfo, std::allocator<RDKit::MolStandardize::ValidationErrorInfo> > {lvalue})
        """
    __instance_size__ = 32
    pass
class Normalizer(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,unsigned int)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, normalizeFilename: str, maxRestarts: int) -> None: ...
    def normalize(self, mol: Mol) -> Mol: 
        """
        normalize( self: Normalizer, mol: Mol) -> Mol

            C++ signature :
                RDKit::ROMol* normalize(RDKit::MolStandardize::Normalizer {lvalue},RDKit::ROMol)
        """
    def normalizeInPlace(self, mol: Mol) -> None: 
        """
        normalizeInPlace( self: Normalizer, mol: Mol) -> None
            modifies the input molecule

            C++ signature :
                void normalizeInPlace(RDKit::MolStandardize::Normalizer {lvalue},RDKit::ROMol {lvalue})
        """
    __instance_size__ = 40
    pass
class RDKitValidation(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    def validate(self, mol: Mol, reportAllFailures: bool = False) -> list: 
        """
        validate( self: RDKitValidation, mol: Mol, reportAllFailures: bool = False) -> list

            C++ signature :
                boost::python::list validate(RDKit::MolStandardize::RDKitValidation {lvalue},RDKit::ROMol [,bool=False])
        """
    __instance_size__ = 32
    pass
class Reionizer(Boost.Python.instance):
    @staticmethod
    @typing.overload
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

            C++ signature :
                void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::vector<RDKit::MolStandardize::ChargeCorrection, std::allocator<RDKit::MolStandardize::ChargeCorrection> >)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str) -> None: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: object, arg2: str, arg3: object) -> None: ...
    def reionize(self, mol: Mol) -> Mol: 
        """
        reionize( self: Reionizer, mol: Mol) -> Mol

            C++ signature :
                RDKit::ROMol* reionize(RDKit::MolStandardize::Reionizer {lvalue},RDKit::ROMol)
        """
    def reionizeInPlace(self, mol: Mol) -> None: 
        """
        reionizeInPlace( self: Reionizer, mol: Mol) -> None
            modifies the input molecule

            C++ signature :
                void reionizeInPlace(RDKit::MolStandardize::Reionizer {lvalue},RDKit::ROMol {lvalue})
        """
    __instance_size__ = 56
    pass
class SmilesTautomerMap(Boost.Python.instance):
    """
    maps SMILES strings to the respective Tautomer objects
    """
    @staticmethod
    def __contains__( arg1: SmilesTautomerMap, arg2: object) -> bool: 
        """
        __contains__( arg1: SmilesTautomerMap, arg2: object) -> bool

            C++ signature :
                bool __contains__(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > > {lvalue},_object*)
        """
    @staticmethod
    def __delitem__( arg1: SmilesTautomerMap, arg2: object) -> None: 
        """
        __delitem__( arg1: SmilesTautomerMap, arg2: object) -> None

            C++ signature :
                void __delitem__(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > > {lvalue},_object*)
        """
    @staticmethod
    def __getitem__( arg1: object, arg2: object) -> object: 
        """
        __getitem__( arg1: object, arg2: object) -> object

            C++ signature :
                boost::python::api::object __getitem__(boost::python::back_reference<std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >&>,_object*)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::_Rb_tree_iterator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > > __iter__(boost::python::back_reference<std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >&>)
        """
    @staticmethod
    def __len__( arg1: SmilesTautomerMap) -> int: 
        """
        __len__( arg1: SmilesTautomerMap) -> int

            C++ signature :
                unsigned long __len__(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > > {lvalue})
        """
    @staticmethod
    def __setitem__( arg1: SmilesTautomerMap, arg2: object, arg3: object) -> None: 
        """
        __setitem__( arg1: SmilesTautomerMap, arg2: object, arg3: object) -> None

            C++ signature :
                void __setitem__(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > > {lvalue},_object*,_object*)
        """
    @staticmethod
    def items( arg1: SmilesTautomerMap) -> tuple: 
        """
        items( arg1: SmilesTautomerMap) -> tuple

            C++ signature :
                boost::python::tuple items(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
    @staticmethod
    def keys( arg1: SmilesTautomerMap) -> tuple: 
        """
        keys( arg1: SmilesTautomerMap) -> tuple

            C++ signature :
                boost::python::tuple keys(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
    @staticmethod
    def values( arg1: SmilesTautomerMap) -> tuple: 
        """
        values( arg1: SmilesTautomerMap) -> tuple

            C++ signature :
                boost::python::tuple values(std::map<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, RDKit::MolStandardize::Tautomer, std::less<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >, std::allocator<std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> > >)
        """
    pass
class Tautomer(Boost.Python.instance):
    """
    used to hold the aromatic and kekulized versions of each tautomer
    """
    @property
    def kekulized(self) -> None:
        """
        kekulized version of the tautomer

        :type: None
        """
    @property
    def tautomer(self) -> None:
        """
        aromatic version of the tautomer

        :type: None
        """
    pass
class TautomerEnumerator(Boost.Python.instance):
    @typing.overload
    def Canonicalize(self, mol: Mol) -> Mol: 
        """
        Canonicalize( self: TautomerEnumerator, mol: Mol) -> Mol
            Returns the canonical tautomer for a molecule.
            
              The default scoring scheme is inspired by the publication:
              M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
              https://doi.org/10.1007/s10822-010-9346-4
            
              Note that the canonical tautomer is very likely not the most stable tautomer
              for any given conditions. The default scoring rules are designed to produce
              "reasonable" tautomers, but the primary concern is that the results are
              canonical: you always get the same canonical tautomer for a molecule
              regardless of what the input tautomer or atom ordering were.

            C++ signature :
                RDKit::ROMol* Canonicalize(RDKit::MolStandardize::TautomerEnumerator,RDKit::ROMol)

            C++ signature :
                RDKit::ROMol* Canonicalize(RDKit::MolStandardize::TautomerEnumerator,RDKit::ROMol,boost::python::api::object)
        """
    @typing.overload
    def Canonicalize(self, mol: Mol, scoreFunc: AtomPairsParameters) -> Mol: ...
    def Enumerate(self, mol: Mol) -> TautomerEnumeratorResult: 
        """
        Enumerate( self: TautomerEnumerator, mol: Mol) -> TautomerEnumeratorResult
            Generates the tautomers for a molecule.
                         
              The enumeration rules are inspired by the publication:
              M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
              https://doi.org/10.1007/s10822-010-9346-4
              
              Note: the definitions used here are that the atoms modified during
              tautomerization are the atoms at the beginning and end of each tautomer
              transform (the H "donor" and H "acceptor" in the transform) and the bonds
              modified during transformation are any bonds whose order is changed during
              the tautomer transform (these are the bonds between the "donor" and the
              "acceptor").

            C++ signature :
                (anonymous namespace)::PyTautomerEnumeratorResult* Enumerate(RDKit::MolStandardize::TautomerEnumerator,RDKit::ROMol)
        """
    @staticmethod
    def GetCallback( arg1: TautomerEnumerator) -> object: 
        """
        GetCallback( arg1: TautomerEnumerator) -> object
            Get the TautomerEnumeratorCallback subclass instance,
            or None if none was set.

            C++ signature :
                boost::python::api::object GetCallback(RDKit::MolStandardize::TautomerEnumerator)
        """
    def GetMaxTautomers(self) -> int: 
        """
        GetMaxTautomers( self: TautomerEnumerator) -> int
            returns the maximum number of tautomers to be generated.

            C++ signature :
                unsigned int GetMaxTautomers(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
    def GetMaxTransforms(self) -> int: 
        """
        GetMaxTransforms( self: TautomerEnumerator) -> int
            returns the maximum number of transformations to be applied.

            C++ signature :
                unsigned int GetMaxTransforms(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
    def GetReassignStereo(self) -> bool: 
        """
        GetReassignStereo( self: TautomerEnumerator) -> bool
            returns whether AssignStereochemistry will be called on each tautomer generated by the Enumerate() method.

            C++ signature :
                bool GetReassignStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
    def GetRemoveBondStereo(self) -> bool: 
        """
        GetRemoveBondStereo( self: TautomerEnumerator) -> bool
            returns whether stereochemistry information will be removed from double bonds involved in tautomerism.

            C++ signature :
                bool GetRemoveBondStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
    def GetRemoveSp3Stereo(self) -> bool: 
        """
        GetRemoveSp3Stereo( self: TautomerEnumerator) -> bool
            returns whether stereochemistry information will be removed from sp3 atoms involved in tautomerism.

            C++ signature :
                bool GetRemoveSp3Stereo(RDKit::MolStandardize::TautomerEnumerator {lvalue})
        """
    @typing.overload
    def PickCanonical(self, typing.Iterable: AtomPairsParameters) -> Mol: 
        """
        PickCanonical( self: TautomerEnumerator, iterable: AtomPairsParameters) -> Mol
            picks the canonical tautomer from an iterable of molecules

            C++ signature :
                RDKit::ROMol* PickCanonical(RDKit::MolStandardize::TautomerEnumerator,boost::python::api::object)

            C++ signature :
                RDKit::ROMol* PickCanonical(RDKit::MolStandardize::TautomerEnumerator,boost::python::api::object,boost::python::api::object)
        """
    @typing.overload
    def PickCanonical(self, typing.Iterable: AtomPairsParameters, scoreFunc: AtomPairsParameters) -> Mol: ...
    @staticmethod
    def ScoreTautomer( mol: Mol) -> int: 
        """
        ScoreTautomer( mol: Mol) -> int
            returns the score for a tautomer using the default scoring scheme.

            C++ signature :
                int ScoreTautomer(RDKit::ROMol)
        """
    @staticmethod
    def SetCallback( arg1: TautomerEnumerator, arg2: object) -> None: 
        """
        SetCallback( arg1: TautomerEnumerator, arg2: object) -> None
            Pass an instance of a class derived from
            TautomerEnumeratorCallback, which must implement the
            __call__() method.

            C++ signature :
                void SetCallback(RDKit::MolStandardize::TautomerEnumerator {lvalue},_object*)
        """
    def SetMaxTautomers(self, maxTautomers: int) -> None: 
        """
        SetMaxTautomers( self: TautomerEnumerator, maxTautomers: int) -> None
            set the maximum number of tautomers to be generated.

            C++ signature :
                void SetMaxTautomers(RDKit::MolStandardize::TautomerEnumerator {lvalue},unsigned int)
        """
    def SetMaxTransforms(self, maxTransforms: int) -> None: 
        """
        SetMaxTransforms( self: TautomerEnumerator, maxTransforms: int) -> None
            set the maximum number of transformations to be applied. This limit is usually hit earlier than the maxTautomers limit and leads to a more linear scaling of CPU time with increasing number of tautomeric centers (see Sitzmann et al.).

            C++ signature :
                void SetMaxTransforms(RDKit::MolStandardize::TautomerEnumerator {lvalue},unsigned int)
        """
    def SetReassignStereo(self, reassignStereo: bool) -> None: 
        """
        SetReassignStereo( self: TautomerEnumerator, reassignStereo: bool) -> None
            set to True if you wish AssignStereochemistry to be called on each tautomer generated by the Enumerate() method. This defaults to True.

            C++ signature :
                void SetReassignStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
    def SetRemoveBondStereo(self, removeBondStereo: bool) -> None: 
        """
        SetRemoveBondStereo( self: TautomerEnumerator, removeBondStereo: bool) -> None
            set to True if you wish stereochemistry information to be removed from double bonds involved in tautomerism. This means that enols will lose their E/Z stereochemistry after going through tautomer enumeration because of the keto-enolic tautomerism. This defaults to True in the RDKit and also in the workflow described by Sitzmann et al.

            C++ signature :
                void SetRemoveBondStereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
    def SetRemoveSp3Stereo(self, removeSp3Stereo: bool) -> None: 
        """
        SetRemoveSp3Stereo( self: TautomerEnumerator, removeSp3Stereo: bool) -> None
            set to True if you wish stereochemistry information to be removed from sp3 atoms involved in tautomerism. This means that S-aminoacids will lose their stereochemistry after going through tautomer enumeration because of the amido-imidol tautomerism. This defaults to True in RDKit, and to False in the workflow described by Sitzmann et al.

            C++ signature :
                void SetRemoveSp3Stereo(RDKit::MolStandardize::TautomerEnumerator {lvalue},bool)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters) -> object: 
        """
        __init__( arg1: AtomPairsParameters) -> object

            C++ signature :
                void* __init__(boost::python::api::object)

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::MolStandardize::CleanupParameters)

            C++ signature :
                void* __init__(boost::python::api::object,RDKit::MolStandardize::TautomerEnumerator)
        """
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, arg2: CleanupParameters) -> object: ...
    @staticmethod
    @typing.overload
    def __init__( arg1: AtomPairsParameters, arg2: TautomerEnumerator) -> object: ...
    tautomerScoreVersion = '1.0.0'
    pass
class TautomerEnumeratorCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
        implement the __call__() method.
        The __call__() method is called in the innermost loop of the
        algorithm, and provides a mechanism to monitor or stop
        its progress.

        To have your callback called, pass an instance of your
        derived class to TautomerEnumerator.SetCallback()
    """
    @staticmethod
    @typing.overload
    def __call__( arg1: TautomerEnumeratorCallback, arg2: Mol, arg3: object) -> bool: 
        """
        __call__( arg1: TautomerEnumeratorCallback, arg2: Mol, arg3: object) -> bool
            This must be implemented in the derived class. Return True if the tautomer enumeration should continue; False if the tautomer enumeration should stop.
            

            C++ signature :
                bool __call__((anonymous namespace)::PyTautomerEnumeratorCallback {lvalue},RDKit::ROMol,RDKit::MolStandardize::TautomerEnumeratorResult)

            C++ signature :
                void __call__((anonymous namespace)::PyTautomerEnumeratorCallback {lvalue},RDKit::ROMol,RDKit::MolStandardize::TautomerEnumeratorResult)
        """
    @staticmethod
    @typing.overload
    def __call__( arg1: TautomerEnumeratorCallback, arg2: Mol, arg3: object) -> None: ...
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    __instance_size__ = 56
    pass
class TautomerEnumeratorResult(Boost.Python.instance):
    """
    used to return tautomer enumeration results
    """
    @staticmethod
    def __call__( arg1: TautomerEnumeratorResult) -> MOL_SPTR_VECT: 
        """
        __call__( arg1: TautomerEnumeratorResult) -> MOL_SPTR_VECT
            tautomers generated by the enumerator

            C++ signature :
                std::vector<boost::shared_ptr<RDKit::ROMol>, std::allocator<boost::shared_ptr<RDKit::ROMol> > > const* __call__((anonymous namespace)::PyTautomerEnumeratorResult {lvalue})
        """
    @staticmethod
    def __getitem__( arg1: TautomerEnumeratorResult, arg2: int) -> Mol: 
        """
        __getitem__( arg1: TautomerEnumeratorResult, arg2: int) -> Mol

            C++ signature :
                RDKit::ROMol* __getitem__((anonymous namespace)::PyTautomerEnumeratorResult {lvalue},int)
        """
    @staticmethod
    def __iter__( arg1: object) -> object: 
        """
        __iter__( arg1: object) -> object

            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, RDKit::MolStandardize::TautomerEnumeratorResult::const_iterator> __iter__(boost::python::back_reference<(anonymous namespace)::PyTautomerEnumeratorResult&>)
        """
    @staticmethod
    def __len__( arg1: TautomerEnumeratorResult) -> int: 
        """
        __len__( arg1: TautomerEnumeratorResult) -> int

            C++ signature :
                int __len__((anonymous namespace)::PyTautomerEnumeratorResult {lvalue})
        """
    @property
    def modifiedAtoms(self) -> None:
        """
        tuple of atom indices modified by the transforms

        :type: None
        """
    @property
    def modifiedBonds(self) -> None:
        """
        tuple of bond indices modified by the transforms

        :type: None
        """
    @property
    def smiles(self) -> None:
        """
        SMILES of tautomers generated by the enumerator

        :type: None
        """
    @property
    def smilesTautomerMap(self) -> None:
        """
        dictionary mapping SMILES strings to the respective Tautomer objects

        :type: None
        """
    @property
    def status(self) -> None:
        """
        whether the enumeration completed or not; see TautomerEnumeratorStatus for possible values

        :type: None
        """
    @property
    def tautomers(self) -> None:
        """
        tautomers generated by the enumerator

        :type: None
        """
    pass
class TautomerEnumeratorStatus(Boost.Python.enum, int):
    Canceled = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled
    Completed = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed
    MaxTautomersReached = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached
    MaxTransformsReached = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached
    __slots__ = ()
    names = {'Completed': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed, 'MaxTautomersReached': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached, 'MaxTransformsReached': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached, 'Canceled': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled}
    values = {0: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed, 1: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached, 2: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached, 3: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled}
    pass
class Uncharger(Boost.Python.instance):
    def __init__(self, canonicalOrder: bool = True) -> None: 
        """
        __init__( self: object, canonicalOrder: bool = True) -> None

            C++ signature :
                void __init__(_object* [,bool=True])
        """
    def uncharge(self, mol: Mol) -> Mol: 
        """
        uncharge( self: Uncharger, mol: Mol) -> Mol

            C++ signature :
                RDKit::ROMol* uncharge(RDKit::MolStandardize::Uncharger {lvalue},RDKit::ROMol)
        """
    def unchargeInPlace(self, mol: Mol) -> None: 
        """
        unchargeInPlace( self: Uncharger, mol: Mol) -> None
            modifies the input molecule

            C++ signature :
                void unchargeInPlace(RDKit::MolStandardize::Uncharger {lvalue},RDKit::ROMol {lvalue})
        """
    __instance_size__ = 96
    pass
class map_indexing_suite_SmilesTautomerMap_entry(Boost.Python.instance):
    @staticmethod
    def __init__( arg1: object) -> None: 
        """
        __init__( arg1: object) -> None

            C++ signature :
                void __init__(_object*)
        """
    @staticmethod
    def __repr__( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> object: 
        """
        __repr__( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> object

            C++ signature :
                boost::python::api::object __repr__(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer>)
        """
    @staticmethod
    def data( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> Tautomer: 
        """
        data( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> Tautomer

            C++ signature :
                RDKit::MolStandardize::Tautomer data(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> {lvalue})
        """
    @staticmethod
    def key( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> str: 
        """
        key( arg1: map_indexing_suite_SmilesTautomerMap_entry) -> str

            C++ signature :
                std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > key(std::pair<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > const, RDKit::MolStandardize::Tautomer> {lvalue})
        """
    __instance_size__ = 112
    pass
def CHARGE_CORRECTIONS() -> object:
    """
    CHARGE_CORRECTIONS() -> object

        C++ signature :
            std::vector<RDKit::MolStandardize::ChargeCorrection, std::allocator<RDKit::MolStandardize::ChargeCorrection> > CHARGE_CORRECTIONS()
    """
def CanonicalTautomer( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    CanonicalTautomer( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Returns the canonical tautomer for the molecule

        C++ signature :
            RDKit::ROMol* CanonicalTautomer(RDKit::ROMol const* [,boost::python::api::object=None])
    """
def ChargeParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    ChargeParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        Returns the uncharged version of the largest fragment

        C++ signature :
            RDKit::ROMol* ChargeParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def Cleanup( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Cleanup( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Standardizes a molecule

        C++ signature :
            RDKit::ROMol* Cleanup(RDKit::ROMol const* [,boost::python::api::object=None])
    """
def CleanupInPlace( mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    CleanupInPlace( mol: Mol, params: AtomPairsParameters = None) -> None
        Standardizes a molecule in place

        C++ signature :
            void CleanupInPlace(RDKit::ROMol* [,boost::python::api::object=None])
    """
def DisconnectOrganometallics( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    DisconnectOrganometallics( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Returns the molecule disconnected using the organometallics rules.

        C++ signature :
            RDKit::ROMol* DisconnectOrganometallics(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
def DisconnectOrganometallicsInPlace( mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    DisconnectOrganometallicsInPlace( mol: Mol, params: AtomPairsParameters = None) -> None
        Disconnects the molecule using the organometallics rules, modifies the input molecule

        C++ signature :
            void DisconnectOrganometallicsInPlace(RDKit::ROMol* [,boost::python::api::object=None])
    """
def FragmentParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    FragmentParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        Returns the largest fragment after doing a cleanup

        C++ signature :
            RDKit::ROMol* FragmentParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def FragmentRemoverFromData( fragmentData: str, leave_last: bool = True, skip_if_all_match: bool = False) -> FragmentRemover:
    """
    FragmentRemoverFromData( fragmentData: str, leave_last: bool = True, skip_if_all_match: bool = False) -> FragmentRemover
        creates a FragmentRemover from a string containing parameter data

        C++ signature :
            RDKit::MolStandardize::FragmentRemover* FragmentRemoverFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=False]])
    """
def GetV1TautomerEnumerator() -> TautomerEnumerator:
    """
    GetV1TautomerEnumerator() -> TautomerEnumerator
        return a TautomerEnumerator using v1 of the enumeration rules

        C++ signature :
            RDKit::MolStandardize::TautomerEnumerator* GetV1TautomerEnumerator()
    """
def IsotopeParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    IsotopeParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        removes all isotopes specifications from the given molecule

        C++ signature :
            RDKit::ROMol* IsotopeParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def Normalize( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Normalize( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Applies a series of standard transformations to correct functional groups and recombine charges

        C++ signature :
            RDKit::ROMol* Normalize(RDKit::ROMol const* [,boost::python::api::object=None])
    """
def NormalizeInPlace( mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    NormalizeInPlace( mol: Mol, params: AtomPairsParameters = None) -> None
        Applies a series of standard transformations to correct functional groups and recombine charges, modifies the input molecule

        C++ signature :
            void NormalizeInPlace(RDKit::ROMol* [,boost::python::api::object=None])
    """
def NormalizerFromData( paramData: str, params: CleanupParameters) -> Normalizer:
    """
    NormalizerFromData( paramData: str, params: CleanupParameters) -> Normalizer
        creates a Normalizer from a string containing normalization SMARTS

        C++ signature :
            RDKit::MolStandardize::Normalizer* NormalizerFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::MolStandardize::CleanupParameters)
    """
def NormalizerFromParams( params: CleanupParameters) -> Normalizer:
    """
    NormalizerFromParams( params: CleanupParameters) -> Normalizer
        creates a Normalizer from CleanupParameters

        C++ signature :
            RDKit::MolStandardize::Normalizer* NormalizerFromParams(RDKit::MolStandardize::CleanupParameters)
    """
def Reionize( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Reionize( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Ensures the strongest acid groups are charged first

        C++ signature :
            RDKit::ROMol* Reionize(RDKit::ROMol const* [,boost::python::api::object=None])
    """
def ReionizeInPlace( mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    ReionizeInPlace( mol: Mol, params: AtomPairsParameters = None) -> None
        Ensures the strongest acid groups are charged first, modifies the input molecule

        C++ signature :
            void ReionizeInPlace(RDKit::ROMol* [,boost::python::api::object=None])
    """
def ReionizerFromData( paramData: str, chargeCorrections: AtomPairsParameters = []) -> Reionizer:
    """
    ReionizerFromData( paramData: str, chargeCorrections: AtomPairsParameters = []) -> Reionizer
        creates a reionizer from a string containing parameter data and a list of charge corrections

        C++ signature :
            RDKit::MolStandardize::Reionizer* ReionizerFromData(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=[]])
    """
def RemoveFragments( mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    RemoveFragments( mol: Mol, params: AtomPairsParameters = None) -> Mol
        Removes fragments from the molecule

        C++ signature :
            RDKit::ROMol* RemoveFragments(RDKit::ROMol const* [,boost::python::api::object=None])
    """
def RemoveFragmentsInPlace( mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    RemoveFragmentsInPlace( mol: Mol, params: AtomPairsParameters = None) -> None
        Removes fragments from the molecule, modifies the input molecule

        C++ signature :
            void RemoveFragmentsInPlace(RDKit::ROMol* [,boost::python::api::object=None])
    """
def StandardizeSmiles( smiles: str) -> str:
    """
    StandardizeSmiles( smiles: str) -> str
        Convenience function for standardizing a SMILES

        C++ signature :
            std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > StandardizeSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def StereoParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    StereoParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        calls removeStereochemistry() on the given molecule

        C++ signature :
            RDKit::ROMol* StereoParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def SuperParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    SuperParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        Returns the super parent. The super parent is the fragment, charge, isotope, stereo, and tautomer parent of the molecule.

        C++ signature :
            RDKit::ROMol* SuperParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def TautomerParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol:
    """
    TautomerParent( mol: Mol, params: AtomPairsParameters = None, skipStandardize: bool = False) -> Mol
        Returns the tautomer parent of a given molecule. The fragment parent is the standardized canonical tautomer of the molecule

        C++ signature :
            RDKit::ROMol* TautomerParent(RDKit::ROMol const* [,boost::python::api::object=None [,bool=False]])
    """
def UpdateParamsFromJSON( arg1: CleanupParameters, arg2: str) -> None:
    """
    UpdateParamsFromJSON( arg1: CleanupParameters, arg2: str) -> None
        updates the cleanup parameters from the provided JSON string

        C++ signature :
            void UpdateParamsFromJSON(RDKit::MolStandardize::CleanupParameters {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
def ValidateSmiles( mol: str) -> list:
    """
    ValidateSmiles( mol: str) -> list

        C++ signature :
            boost::python::list ValidateSmiles(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
