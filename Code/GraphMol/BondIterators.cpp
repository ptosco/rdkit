// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "BondIterators.h"

namespace RDKit {

BondIterator_::BondIterator_(ROMol *mol) {
  _mol = mol;
  boost::tie(_beg, _end) = mol->getEdges();
  _pos = _beg;
  std::cerr << "1) BondIterator_::BondIterator_ this " << this << ", _mol " << _mol << std::endl;
};
BondIterator_::BondIterator_(ROMol *mol, ROMol::EDGE_ITER pos) {
  _mol = mol;
  boost::tie(_beg, _end) = mol->getEdges();

  _pos = pos;
  std::cerr << "2) BondIterator_::BondIterator_ this " << this << ", _mol " << _mol << std::endl;
};
BondIterator_::BondIterator_(const BondIterator_ &other) {
  _mol = other._mol;
  _pos = other._pos;
  _beg = other._beg;
  _end = other._end;
  std::cerr << "3) BondIterator_::BondIterator_ this " << this << ", _mol " << _mol << std::endl;
}

BondIterator_ &BondIterator_::operator=(const BondIterator_ &other) {
  _mol = other._mol;
  _pos = other._pos;
  _beg = other._beg;
  _end = other._end;
  return *this;
}

bool BondIterator_::operator==(const BondIterator_ &other) const {
  return _mol == other._mol && _pos == other._pos;
}
bool BondIterator_::operator!=(const BondIterator_ &other) const {
  return _mol != other._mol || _pos != other._pos;
}

Bond *BondIterator_::operator*() const {
  std::cerr << "1) BondIterator_::operator* this " << this << ", _mol " << _mol << ", _pos == _end " << (_pos == _end) << std::endl;
  auto res = (*_mol)[*_pos];
  std::cerr << "2) BondIterator_::operator* this " << this << ", _mol " << _mol << ", _pos == _end " << (_pos == _end) << std::endl;
  return res;
}
// pre-increment
BondIterator_ &BondIterator_::operator++() {
  PRECONDITION(_pos != _end, "bad initial position")
  _pos++;
  return *this;
}
BondIterator_ BondIterator_::operator++(int) {
  PRECONDITION(_pos != _end, "bad initial position")
  BondIterator_ res(*this);
  _pos++;
  return res;
}
// pre-decrement
BondIterator_ &BondIterator_::operator--() {
  if (_pos == _beg) {
    _pos = _end;
  } else {
    _pos--;
  }
  return *this;
}
BondIterator_ BondIterator_::operator--(int) {
  BondIterator_ res(*this);
  if (_pos == _beg) {
    _pos = _end;
  } else {
    _pos--;
  }
  return res;
}

// CONST
ConstBondIterator_::ConstBondIterator_(ROMol const *mol) {
  _mol = mol;
  boost::tie(_beg, _end) = mol->getEdges();
  _pos = _beg;
};
ConstBondIterator_::ConstBondIterator_(ROMol const *mol, ROMol::EDGE_ITER pos) {
  _mol = mol;
  boost::tie(_beg, _end) = mol->getEdges();

  _pos = pos;
};
ConstBondIterator_::ConstBondIterator_(const ConstBondIterator_ &other) {
  _mol = other._mol;
  _pos = other._pos;
  _beg = other._beg;
  _end = other._end;
}
ConstBondIterator_ &ConstBondIterator_::operator=(
    const ConstBondIterator_ &other) {
  _mol = other._mol;
  _pos = other._pos;
  _beg = other._beg;
  _end = other._end;
  return *this;
}
bool ConstBondIterator_::operator==(const ConstBondIterator_ &other) const {
  return _mol == other._mol && _pos == other._pos;
}
bool ConstBondIterator_::operator!=(const ConstBondIterator_ &other) const {
  return _mol != other._mol || _pos != other._pos;
}

Bond const *ConstBondIterator_::operator*() const { return (*_mol)[*_pos]; }
// pre-increment
ConstBondIterator_ &ConstBondIterator_::operator++() {
  PRECONDITION(_pos != _end, "bad initial position")
  _pos++;
  return *this;
}
ConstBondIterator_ ConstBondIterator_::operator++(int) {
  ConstBondIterator_ res(*this);
  PRECONDITION(_pos != _end, "bad initial position")
  _pos++;
  return res;
}
// pre-decrement
ConstBondIterator_ &ConstBondIterator_::operator--() {
  if (_pos == _beg) {
    _pos = _end;
  } else {
    _pos--;
  }
  return *this;
}
ConstBondIterator_ ConstBondIterator_::operator--(int) {
  ConstBondIterator_ res(*this);
  if (_pos == _beg) {
    _pos = _end;
  } else {
    _pos--;
  }
  return res;
}
}  // namespace RDKit
