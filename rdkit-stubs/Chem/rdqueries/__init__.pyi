"""Module containing RDKit functionality for querying molecules."""
from __future__ import annotations
import rdkit.Chem.rdqueries
import typing

__all__ = [
    "AAtomQueryAtom",
    "AHAtomQueryAtom",
    "AtomNumEqualsQueryAtom",
    "AtomNumGreaterQueryAtom",
    "AtomNumLessQueryAtom",
    "ExplicitDegreeEqualsQueryAtom",
    "ExplicitDegreeGreaterQueryAtom",
    "ExplicitDegreeLessQueryAtom",
    "ExplicitValenceEqualsQueryAtom",
    "ExplicitValenceGreaterQueryAtom",
    "ExplicitValenceLessQueryAtom",
    "FormalChargeEqualsQueryAtom",
    "FormalChargeGreaterQueryAtom",
    "FormalChargeLessQueryAtom",
    "HCountEqualsQueryAtom",
    "HCountGreaterQueryAtom",
    "HCountLessQueryAtom",
    "HasBitVectPropWithValueQueryAtom",
    "HasBoolPropWithValueQueryAtom",
    "HasBoolPropWithValueQueryBond",
    "HasChiralTagQueryAtom",
    "HasDoublePropWithValueQueryAtom",
    "HasDoublePropWithValueQueryBond",
    "HasIntPropWithValueQueryAtom",
    "HasIntPropWithValueQueryBond",
    "HasPropQueryAtom",
    "HasPropQueryBond",
    "HasStringPropWithValueQueryAtom",
    "HasStringPropWithValueQueryBond",
    "HybridizationEqualsQueryAtom",
    "HybridizationGreaterQueryAtom",
    "HybridizationLessQueryAtom",
    "InNRingsEqualsQueryAtom",
    "InNRingsGreaterQueryAtom",
    "InNRingsLessQueryAtom",
    "IsAliphaticQueryAtom",
    "IsAromaticQueryAtom",
    "IsBridgeheadQueryAtom",
    "IsInRingQueryAtom",
    "IsUnsaturatedQueryAtom",
    "IsotopeEqualsQueryAtom",
    "IsotopeGreaterQueryAtom",
    "IsotopeLessQueryAtom",
    "MAtomQueryAtom",
    "MHAtomQueryAtom",
    "MassEqualsQueryAtom",
    "MassGreaterQueryAtom",
    "MassLessQueryAtom",
    "MinRingSizeEqualsQueryAtom",
    "MinRingSizeGreaterQueryAtom",
    "MinRingSizeLessQueryAtom",
    "MissingChiralTagQueryAtom",
    "NonHydrogenDegreeEqualsQueryAtom",
    "NonHydrogenDegreeGreaterQueryAtom",
    "NonHydrogenDegreeLessQueryAtom",
    "NumAliphaticHeteroatomNeighborsEqualsQueryAtom",
    "NumAliphaticHeteroatomNeighborsGreaterQueryAtom",
    "NumAliphaticHeteroatomNeighborsLessQueryAtom",
    "NumHeteroatomNeighborsEqualsQueryAtom",
    "NumHeteroatomNeighborsGreaterQueryAtom",
    "NumHeteroatomNeighborsLessQueryAtom",
    "NumRadicalElectronsEqualsQueryAtom",
    "NumRadicalElectronsGreaterQueryAtom",
    "NumRadicalElectronsLessQueryAtom",
    "QAtomQueryAtom",
    "QHAtomQueryAtom",
    "RingBondCountEqualsQueryAtom",
    "RingBondCountGreaterQueryAtom",
    "RingBondCountLessQueryAtom",
    "TotalDegreeEqualsQueryAtom",
    "TotalDegreeGreaterQueryAtom",
    "TotalDegreeLessQueryAtom",
    "TotalValenceEqualsQueryAtom",
    "TotalValenceGreaterQueryAtom",
    "TotalValenceLessQueryAtom",
    "XAtomQueryAtom",
    "XHAtomQueryAtom"
]


def AAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    AAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when AAtom is True.

        C++ signature :
            RDKit::QueryAtom* AAtomQueryAtom([ bool=False])
    """
def AHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    AHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when AHAtom is True.

        C++ signature :
            RDKit::QueryAtom* AHAtomQueryAtom([ bool=False])
    """
def AtomNumEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* AtomNumEqualsQueryAtom(int [,bool=False])
    """
def AtomNumGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* AtomNumGreaterQueryAtom(int [,bool=False])
    """
def AtomNumLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    AtomNumLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where AtomNum is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* AtomNumLessQueryAtom(int [,bool=False])
    """
def ExplicitDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeEqualsQueryAtom(int [,bool=False])
    """
def ExplicitDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeGreaterQueryAtom(int [,bool=False])
    """
def ExplicitDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitDegreeLessQueryAtom(int [,bool=False])
    """
def ExplicitValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceEqualsQueryAtom(int [,bool=False])
    """
def ExplicitValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceGreaterQueryAtom(int [,bool=False])
    """
def ExplicitValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    ExplicitValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where ExplicitValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* ExplicitValenceLessQueryAtom(int [,bool=False])
    """
def FormalChargeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* FormalChargeEqualsQueryAtom(int [,bool=False])
    """
def FormalChargeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* FormalChargeGreaterQueryAtom(int [,bool=False])
    """
def FormalChargeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    FormalChargeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where FormalCharge is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* FormalChargeLessQueryAtom(int [,bool=False])
    """
def HCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* HCountEqualsQueryAtom(int [,bool=False])
    """
def HCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HCountGreaterQueryAtom(int [,bool=False])
    """
def HCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where HCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HCountLessQueryAtom(int [,bool=False])
    """
def HasBitVectPropWithValueQueryAtom( propname: str, val: ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> QueryAtom:
    """
    HasBitVectPropWithValueQueryAtom( propname: str, val: ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified explicit bit vector value.  The Tolerance is the allowed Tanimoto difference

        C++ signature :
            RDKit::QueryAtom* HasBitVectPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,ExplicitBitVect [,bool=False [,float=0]])
    """
def HasBoolPropWithValueQueryAtom( propname: str, val: bool, negate: bool = False) -> QueryAtom:
    """
    HasBoolPropWithValueQueryAtom( propname: str, val: bool, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified boolean value.

        C++ signature :
            RDKit::QueryAtom* HasBoolPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
def HasBoolPropWithValueQueryBond( propname: str, val: bool, negate: bool = False) -> QueryBond:
    """
    HasBoolPropWithValueQueryBond( propname: str, val: bool, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified boolean value.

        C++ signature :
            RDKit::QueryBond* HasBoolPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,bool [,bool=False])
    """
def HasChiralTagQueryAtom(negate: bool = False) -> QueryAtom:
    """
    HasChiralTagQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when HasChiralTag is True.

        C++ signature :
            RDKit::QueryAtom* HasChiralTagQueryAtom([ bool=False])
    """
def HasDoublePropWithValueQueryAtom( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryAtom:
    """
    HasDoublePropWithValueQueryAtom( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified value +- tolerance

        C++ signature :
            RDKit::QueryAtom* HasDoublePropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
def HasDoublePropWithValueQueryBond( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryBond:
    """
    HasDoublePropWithValueQueryBond( propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified value +- tolerance

        C++ signature :
            RDKit::QueryBond* HasDoublePropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,double [,bool=False [,double=0.0]])
    """
def HasIntPropWithValueQueryAtom( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryAtom:
    """
    HasIntPropWithValueQueryAtom( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified int value.

        C++ signature :
            RDKit::QueryAtom* HasIntPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
def HasIntPropWithValueQueryBond( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryBond:
    """
    HasIntPropWithValueQueryBond( propname: str, val: int, negate: bool = False, tolerance: int = 0) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified int value.

        C++ signature :
            RDKit::QueryBond* HasIntPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,int [,bool=False [,int=0]])
    """
def HasPropQueryAtom( propname: str, negate: bool = False) -> QueryAtom:
    """
    HasPropQueryAtom( propname: str, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' exists in the atom.

        C++ signature :
            RDKit::QueryAtom* HasPropQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasPropQueryBond( propname: str, negate: bool = False) -> QueryBond:
    """
    HasPropQueryBond( propname: str, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' exists in the bond.

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])

        C++ signature :
            RDKit::QueryBond* HasPropQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasStringPropWithValueQueryAtom( propname: str, val: str, negate: bool = False) -> QueryAtom:
    """
    HasStringPropWithValueQueryAtom( propname: str, val: str, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches when the propery 'propname' has the specified string value.

        C++ signature :
            RDKit::QueryAtom* HasStringPropWithValueQueryAtom(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HasStringPropWithValueQueryBond( propname: str, val: str, negate: bool = False) -> QueryBond:
    """
    HasStringPropWithValueQueryBond( propname: str, val: str, negate: bool = False) -> QueryBond
        Returns a QueryBond that matches when the propery 'propname' has the specified string value.

        C++ signature :
            RDKit::QueryBond* HasStringPropWithValueQueryBond(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False])
    """
def HybridizationEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* HybridizationEqualsQueryAtom(int [,bool=False])
    """
def HybridizationGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HybridizationGreaterQueryAtom(int [,bool=False])
    """
def HybridizationLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    HybridizationLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Hybridization is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* HybridizationLessQueryAtom(int [,bool=False])
    """
def InNRingsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* InNRingsEqualsQueryAtom(int [,bool=False])
    """
def InNRingsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* InNRingsGreaterQueryAtom(int [,bool=False])
    """
def InNRingsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    InNRingsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where InNRings is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* InNRingsLessQueryAtom(int [,bool=False])
    """
def IsAliphaticQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsAliphaticQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsAliphatic is True.

        C++ signature :
            RDKit::QueryAtom* IsAliphaticQueryAtom([ bool=False])
    """
def IsAromaticQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsAromaticQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsAromatic is True.

        C++ signature :
            RDKit::QueryAtom* IsAromaticQueryAtom([ bool=False])
    """
def IsBridgeheadQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsBridgeheadQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsBridgehead is True.

        C++ signature :
            RDKit::QueryAtom* IsBridgeheadQueryAtom([ bool=False])
    """
def IsInRingQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsInRingQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsInRing is True.

        C++ signature :
            RDKit::QueryAtom* IsInRingQueryAtom([ bool=False])
    """
def IsUnsaturatedQueryAtom(negate: bool = False) -> QueryAtom:
    """
    IsUnsaturatedQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when IsUnsaturated is True.

        C++ signature :
            RDKit::QueryAtom* IsUnsaturatedQueryAtom([ bool=False])
    """
def IsotopeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* IsotopeEqualsQueryAtom(int [,bool=False])
    """
def IsotopeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* IsotopeGreaterQueryAtom(int [,bool=False])
    """
def IsotopeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    IsotopeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Isotope is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* IsotopeLessQueryAtom(int [,bool=False])
    """
def MAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MAtom is True.

        C++ signature :
            RDKit::QueryAtom* MAtomQueryAtom([ bool=False])
    """
def MHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MHAtom is True.

        C++ signature :
            RDKit::QueryAtom* MHAtomQueryAtom([ bool=False])
    """
def MassEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* MassEqualsQueryAtom(int [,bool=False])
    """
def MassGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MassGreaterQueryAtom(int [,bool=False])
    """
def MassLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MassLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where Mass is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MassLessQueryAtom(int [,bool=False])
    """
def MinRingSizeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* MinRingSizeEqualsQueryAtom(int [,bool=False])
    """
def MinRingSizeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MinRingSizeGreaterQueryAtom(int [,bool=False])
    """
def MinRingSizeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    MinRingSizeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where MinRingSize is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* MinRingSizeLessQueryAtom(int [,bool=False])
    """
def MissingChiralTagQueryAtom(negate: bool = False) -> QueryAtom:
    """
    MissingChiralTagQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when MissingChiralTag is True.

        C++ signature :
            RDKit::QueryAtom* MissingChiralTagQueryAtom([ bool=False])
    """
def NonHydrogenDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeEqualsQueryAtom(int [,bool=False])
    """
def NonHydrogenDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeGreaterQueryAtom(int [,bool=False])
    """
def NonHydrogenDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NonHydrogenDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NonHydrogenDegreeLessQueryAtom(int [,bool=False])
    """
def NumAliphaticHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
def NumAliphaticHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
def NumAliphaticHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumAliphaticHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumAliphaticHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsEqualsQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsGreaterQueryAtom(int [,bool=False])
    """
def NumHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumHeteroatomNeighborsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumHeteroatomNeighborsLessQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsEqualsQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsGreaterQueryAtom(int [,bool=False])
    """
def NumRadicalElectronsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    NumRadicalElectronsLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* NumRadicalElectronsLessQueryAtom(int [,bool=False])
    """
def QAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    QAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when QAtom is True.

        C++ signature :
            RDKit::QueryAtom* QAtomQueryAtom([ bool=False])
    """
def QHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    QHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when QHAtom is True.

        C++ signature :
            RDKit::QueryAtom* QHAtomQueryAtom([ bool=False])
    """
def RingBondCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* RingBondCountEqualsQueryAtom(int [,bool=False])
    """
def RingBondCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* RingBondCountGreaterQueryAtom(int [,bool=False])
    """
def RingBondCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    RingBondCountLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where RingBondCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* RingBondCountLessQueryAtom(int [,bool=False])
    """
def TotalDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* TotalDegreeEqualsQueryAtom(int [,bool=False])
    """
def TotalDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalDegreeGreaterQueryAtom(int [,bool=False])
    """
def TotalDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalDegreeLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalDegreeLessQueryAtom(int [,bool=False])
    """
def TotalValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceEqualsQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.

        C++ signature :
            RDKit::QueryAtom* TotalValenceEqualsQueryAtom(int [,bool=False])
    """
def TotalValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceGreaterQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalValenceGreaterQueryAtom(int [,bool=False])
    """
def TotalValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom:
    """
    TotalValenceLessQueryAtom( val: int, negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms where TotalValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API

        C++ signature :
            RDKit::QueryAtom* TotalValenceLessQueryAtom(int [,bool=False])
    """
def XAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    XAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when XAtom is True.

        C++ signature :
            RDKit::QueryAtom* XAtomQueryAtom([ bool=False])
    """
def XHAtomQueryAtom(negate: bool = False) -> QueryAtom:
    """
    XHAtomQueryAtom(negate: bool = False) -> QueryAtom
        Returns a QueryAtom that matches atoms when XHAtom is True.

        C++ signature :
            RDKit::QueryAtom* XHAtomQueryAtom([ bool=False])
    """
