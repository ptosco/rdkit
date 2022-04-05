//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#include <list>
#include <algorithm>
#include <cmath>
#include "SubstructMatchCustom.h"

namespace RDKit {
namespace FMCS {
class RDKIT_FMCS_EXPORT RingMatchTableSet {
  class RingMatchTable {
    FMCS::MatchTable MatchMatrix;
    std::map<const INT_VECT*, unsigned int> RingIndex;

   public:
    inline void clear() {
      MatchMatrix.clear();
      RingIndex.clear();
    }
    inline void resize(unsigned int s1, unsigned int s2) {
      MatchMatrix.resize(s1, s2);
      for (size_t i = 0; i < s1; i++) {
        for (size_t j = 0; j < s2; j++) {
          MatchMatrix.set(i, j, false);
        }
      }
    }
    inline void makeRingIndex(const ROMol* mol2) {
      unsigned int i = 0;
      // for each TARGET ring
      const auto& rings2 = mol2->getRingInfo()->bondRings();
      for (const auto &r2 : rings2) {
        RingIndex[&r2] = i++;
      }
    }
    inline bool isEqual(unsigned int i, const INT_VECT* r2) const {
      auto ri = getRingIndex(r2);
      return (ri == -1 ? false : MatchMatrix.at(i, ri));
    }
    inline void setMatch(unsigned int i, const INT_VECT* r2) {
      auto ri = getRingIndex(r2);
      CHECK_INVARIANT(ri != -1, "Failed to find ring in RingIndex");
      MatchMatrix.set(i, ri, true);
    }

   private:
    inline int getRingIndex(const INT_VECT* r2) const {
      auto j = RingIndex.find(r2);
      return (j == RingIndex.end() ? -1 : j->second);
    }
  };

 private:
  std::vector<std::vector<size_t>>* QueryBondRingsIndeces{nullptr};
  std::map<const ROMol*, std::vector<std::vector<size_t>>>
      TargetBondRingsIndecesSet;  // by target molecules

  std::map<const ROMol*, RingMatchTable> MatchMatrixSet;  // by target molecules
  std::map<const INT_VECT*, unsigned int> QueryRingIndex;

 public:
  RingMatchTableSet() {}

  inline void clear() {
    if (QueryBondRingsIndeces) {
      QueryBondRingsIndeces->clear();
    }
    TargetBondRingsIndecesSet.clear();
    MatchMatrixSet.clear();
    QueryRingIndex.clear();
  }

  inline bool isQueryBondInRing(unsigned int bi) const {
    return getQueryBondRings(bi).empty();
  }
  inline const std::vector<size_t>& getQueryBondRings(unsigned int bi) const {
    return (*QueryBondRingsIndeces)[bi];
  }

  inline bool isTargetBondInRing(const ROMol* target, unsigned int bi) const {
    return getTargetBondRings(target, bi).empty();
  }

  inline const std::vector<size_t>& getTargetBondRings(const ROMol* target,
                                                       unsigned int bi) const {
    auto i = TargetBondRingsIndecesSet.find(target);
    CHECK_INVARIANT(i != TargetBondRingsIndecesSet.end(), "Failed to find target in TargetBondRingsIndecesSet");
    return i->second[bi];
  }

  inline bool isEqual(const INT_VECT* r1, const INT_VECT* r2,
                      const ROMol* mol2) const {
    const RingMatchTable& m = getTargetMatchMatrix(mol2);
    auto i = getQueryRingIndex(r1);
    return (i == -1 ? false : m.isEqual(i, r2));
  }

  void init(const ROMol* query) {
    MatchMatrixSet.clear();
    // fill out QueryRingIndex
    unsigned int i = 0;
    const auto& rings = query->getRingInfo()->bondRings();
    for (const auto &r : rings) {
      QueryRingIndex[&r] = i++;
    }
    TargetBondRingsIndecesSet.clear();
    QueryBondRingsIndeces = &TargetBondRingsIndecesSet[query];
    QueryBondRingsIndeces->resize(query->getNumBonds());
    size_t ri = 0;
    for (auto r = rings.begin(); r != rings.end(); r++, ri++) {
      for (auto bi : *r) {  // all bonds in the ring
        (*QueryBondRingsIndeces)[bi].push_back(ri);
      }
    }
  }
  inline void addTargetBondRingsIndeces(const ROMol* mol2) {
    auto &m = TargetBondRingsIndecesSet[mol2];
    m.resize(mol2->getNumBonds());

    size_t ri = 0;
    const auto &rings = mol2->getRingInfo()->bondRings();
    for (auto r = rings.begin(); r != rings.end(); r++, ri++) {
      for (auto bi : *r) {  // all bonds in the ring
        m[bi].push_back(ri);
      }
    }
  }

  void computeRingMatchTable(
      const ROMol* query, const ROMol* targetMolecule,
      const MCSParameters& parameters) {  // call it for all targets
    const auto& rings1 = query->getRingInfo()->bondRings();
    const auto& rings2 = targetMolecule->getRingInfo()->bondRings();
    auto &m = addTargetMatchMatrix(targetMolecule, rings1.size(), rings2.size());
    unsigned int i = 0;
    // for each query ring
    for (auto r1 = rings1.begin(); r1 != rings1.end(); r1++, i++) {
      FMCS::Graph graph1;
      makeRingGraph(graph1, *r1, query);  // for each query ring bond ADD all atoms and bonds

      // for each TARGET ring
      for (auto r2 : rings2) {
        if (r1->size() != r2.size()) {  // rings are different
          continue;
        }
        FMCS::Graph graph2;
        // for each TAG ring bond ADD all atoms and bonds
        makeRingGraph(graph2, r2, targetMolecule);

        // check ring substruct match
        MCSBondCompareParameters bp = parameters.BondCompareParameters;
        bp.RingMatchesRingOnly = false;
        bp.CompleteRingsOnly = false;
        bool match =
#ifdef NEVER_xxx_PRECOMPUTED_TABLES_MATCH  // not computed yet, because
                                           // MatchTable computation uses this
                                           // ring info table
            FMCS::SubstructMatchCustomTable(graph2, graph1, tag->AtomMatchTable,
                                            tag->BondMatchTable);
#else  // noticeable slowly:
            FMCS::SubstructMatchCustom(
                graph2, *targetMolecule, graph1, *query, parameters.AtomTyper,
                parameters.BondTyper, nullptr, parameters.AtomCompareParameters,
                bp, parameters.CompareFunctionsUserData);
#endif
        if (match) {
          m.setMatch(i, &r2);
        }
      }
    }
  }

 private:
  void makeRingGraph(FMCS::Graph& g, const INT_VECT& ring,
                     const ROMol* mol) const {  // ADD all atoms and bonds
    std::map<const Atom*, unsigned int> atomMap;

    for (size_t i = 0; i < ring.size(); i++) {
      const auto bond = mol->getBondWithIdx(ring[i]);
      const auto atom1 = bond->getBeginAtom();
      const auto atom2 = bond->getEndAtom();
      unsigned int j1 = NotSet;
      unsigned int j2 = NotSet;
      auto ai = atomMap.find(atom1);
      if (atomMap.end() != ai) {
        j1 = ai->second;
      }
      ai = atomMap.find(atom2);
      if (atomMap.end() != ai) {
        j2 = ai->second;
      }
      if (NotSet == j1) {
        j1 = g.m_vertices.size();
        atomMap[atom1] = j1;
        g.addAtom(atom1->getIdx());
      }
      if (NotSet == j2) {
        j2 = g.m_vertices.size();
        atomMap[atom2] = j2;
        g.addAtom(atom2->getIdx());
      }
      g.addBond(ring[i], j1, j2);
    }
  }

  inline int getQueryRingIndex(const INT_VECT* r1) const {
    auto i = QueryRingIndex.find(r1);
    return (i == QueryRingIndex.end() ? -1 : i->second);
  }
  inline const RingMatchTable& getTargetMatchMatrix(const ROMol* mol2) const {
    auto mi = MatchMatrixSet.find(mol2);
    CHECK_INVARIANT(mi != MatchMatrixSet.end(), "Failed to find mol2 in MatchMatrixSet");
    return mi->second;
  }

  inline RingMatchTable& addTargetMatchMatrix(const ROMol* mol2, unsigned int s1,
                                              unsigned int s2) {
    auto &m = MatchMatrixSet[mol2];
    m.resize(s1, s2);
    m.makeRingIndex(mol2);
    return m;
  }
};
}  // namespace FMCS
}  // namespace RDKit
