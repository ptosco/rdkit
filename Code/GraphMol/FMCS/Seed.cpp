//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MaximumCommonSubgraph.h"
#include "Composition2N.h"
#include "Seed.h"

#include "DebugTrace.h"
#include "../SmilesParse/SmilesWrite.h"

#include <set>

//#define DEBUG_FMCS 1

namespace RDKit {
namespace FMCS {

unsigned Seed::addAtom(const Atom* atom) {
  unsigned i = MoleculeFragment.AtomsIdx.size();
  unsigned aqi = atom->getIdx();
  MoleculeFragment.Atoms.push_back(atom);
  MoleculeFragment.AtomsIdx.push_back(aqi);
  MoleculeFragment.SeedAtomIdxMap[aqi] = i;
  Topology.addAtom(aqi);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addAtom(aqi);
#endif
  return i;
}

unsigned Seed::addBond(const Bond* bond) {
  unsigned b = bond->getIdx();
  if (ExcludedBonds[b]) {
    throw -1;  // never, check the implementation
  }
  ExcludedBonds[b] = true;
  MoleculeFragment.BondsIdx.push_back(b);
  MoleculeFragment.Bonds.push_back(bond);
  // remap idx to seed's indices:
  unsigned i = MoleculeFragment.SeedAtomIdxMap[bond->getBeginAtomIdx()];
  unsigned j = MoleculeFragment.SeedAtomIdxMap[bond->getEndAtomIdx()];
  Topology.addBond(b, i, j);
#ifdef DUP_SUBSTRUCT_CACHE
  DupCacheKey.addBond(b);
#endif
  return getNumBonds();
}

void Seed::fillNewBonds(const ROMol& qmol) {
  auto excludedBonds = ExcludedBonds;
  for (unsigned srcAtomIdx = LastAddedAtomsBeginIdx; srcAtomIdx < getNumAtoms();
       srcAtomIdx++) {  // all atoms added on previous growing only
    const Atom* atom = MoleculeFragment.Atoms[srcAtomIdx];
    ROMol::OEDGE_ITER beg, end;
    for (boost::tie(beg, end) = qmol.getAtomBonds(atom); beg != end;
         beg++) {  // all bonds from MoleculeFragment.Atoms[srcAtomIdx]
      const Bond* bond = &*(qmol[*beg]);
      if (!excludedBonds.test(bond->getIdx())) {
        // already in the seed or NewBonds list from another atom in a RING
        excludedBonds.set(bond->getIdx());
        unsigned ai = (atom == bond->getBeginAtom()) ? bond->getEndAtomIdx()
                                                     : bond->getBeginAtomIdx();
        const Atom* end_atom = qmol.getAtomWithIdx(ai);
        unsigned end_atom_idx = NotSet;
        for (unsigned i = 0; i < getNumAtoms(); i++) {
          if (end_atom ==
              MoleculeFragment.Atoms[i]) {  // already exists in this seed
            end_atom_idx = i;
            break;
          }
        }
        NewBonds.emplace_back(srcAtomIdx, bond->getIdx(), ai, end_atom_idx,
                              NotSet == end_atom_idx ? end_atom : nullptr);
      }
    }
  }
}

void Seed::grow(MaximumCommonSubgraph& mcs) const {
  const ROMol& qmol = mcs.getQueryMolecule();
  std::set<unsigned int> newAtomsSet;  // keep track of newly added atoms

  if (!canGrowBiggerThan(mcs.getMaxNumberBonds(),
                         mcs.getMaxNumberAtoms())) {  // prune() parent
    GrowingStage = NotSet;                            // finished
#ifdef VERBOSE_STATISTICS_ON
    ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
    return;
  }

  if (0 == GrowingStage) {
    // 0. Fill out list of all directly connected outgoing bonds
    ((Seed*)this)
        ->fillNewBonds(
            qmol);  // non const method, multistage growing optimisation

    if (NewBonds.empty()) {
      GrowingStage = NotSet;  // finished
      return;
    }

    // 1. Check and add the biggest child seed with all outgoing bonds added:
    // Add all bonds at first (build the biggest child seed). All new atoms are
    // already in the seed
    Seed seed;
    seed.createFromParent(this);

#ifdef DEBUG_FMCS
    std::cerr << "Seed::grow NewBonds.size() " << NewBonds.size() << std::endl;
#endif
    for (std::vector<NewBond>::const_iterator nbi = NewBonds.begin();
         nbi != NewBonds.end(); nbi++) {
      unsigned aIdx = nbi->EndAtomIdx;
      if (NotSet == aIdx) {  // new atom
        // check if new bonds simultaneously close a ring
        if (newAtomsSet.find(nbi->NewAtomIdx) == newAtomsSet.end()) {
          const Atom* end_atom = nbi->NewAtom;
          aIdx = seed.addAtom(end_atom);
          newAtomsSet.insert(nbi->NewAtomIdx);
        }
      }
      const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
      seed.addBond(src_bond);
    }
#ifdef VERBOSE_STATISTICS_ON
    { ++mcs.VerboseStatistics.Seed; }
#endif
    seed.RemainingBonds = RemainingBonds - NewBonds.size();  // Added ALL !!!
    seed.RemainingAtoms =
        RemainingAtoms - newAtomsSet.size();  // new atoms added to seed

#ifdef DEBUG_FMCS
    std::cerr << "Seed::grow seed.RemainingBonds " << seed.RemainingBonds << " seed.RemainingAtoms " << seed.RemainingAtoms << std::endl;
#endif
    // prune() Best Sizes
    if (!seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                                mcs.getMaxNumberAtoms())) {
      GrowingStage = NotSet;
#ifdef VERBOSE_STATISTICS_ON
      ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
      return;  // the biggest possible subgraph from this seed is too small for
               // future growing. So, skip ALL children !
    }
    seed.MatchResult = MatchResult;
    bool allMatched = mcs.checkIfMatchAndAppend(
        seed);  // this seed + all extern bonds is a part of MCS
#ifdef DEBUG_FMCS
    std::cerr << "Seed::grow allMatched " << allMatched << std::endl;
#endif

    GrowingStage = 1;
    if (allMatched && NewBonds.size() > 1) {
#ifdef DEBUG_FMCS
      std::cerr << "1) Seed::grow return" << std::endl;
#endif
      return;  // grow deep first. postpone next growing steps
    }
  }
  // 2. Check and add all 2^N-1-1 other possible seeds:
  if (1 == NewBonds.size()) {
    GrowingStage = NotSet;
#ifdef DEBUG_FMCS
    std::cerr << "2) Seed::grow return" << std::endl;
#endif
    return;  // everything has been done
  }
  // OPTIMISATION:
  // check each single bond first: if (this seed + single bond) does not exist
  // in MCS, exclude this new bond from growing this seed.
  unsigned numErasedNewBonds = 0;
  for (auto& nbi : NewBonds) {
#ifdef VERBOSE_STATISTICS_ON
    { ++mcs.VerboseStatistics.Seed; }
#endif
    Seed seed;
    seed.createFromParent(this);

    // existed in this parent seed (ring) or -1
    unsigned aIdx = nbi.EndAtomIdx;
    if (NotSet == aIdx) {  // new atom
      const Atom* end_atom = nbi.NewAtom;
      aIdx = seed.addAtom(end_atom);
    }
    const Bond* src_bond = qmol.getBondWithIdx(nbi.BondIdx);
    seed.addBond(src_bond);
    seed.computeRemainingSize(qmol);

    if (seed.canGrowBiggerThan(mcs.getMaxNumberBonds(),
                               mcs.getMaxNumberAtoms())) {  // prune()
      if (!MatchResult.empty()) {
#ifdef DEBUG_FMCS
        std::cerr << "3) Seed::grow here" << std::endl;
#endif
        seed.MatchResult = MatchResult;
      }
      if (!mcs.checkIfMatchAndAppend(seed)) {
#ifdef DEBUG_FMCS
        std::cerr << "4) Seed::grow here" << std::endl;
#endif
        nbi.BondIdx = NotSet;  // exclude this new bond from growing this seed
                               // - decrease 2^^N-1 to 2^^k-1, k<N.
        ++numErasedNewBonds;
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.SingleBondExcluded;
#endif
      }
    } else {  // seed too small
#ifdef VERBOSE_STATISTICS_ON
      ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
    }
  }
#ifdef DEBUG_FMCS
  std::cerr << "5) Seed::grow here" << std::endl;
#endif

  if (numErasedNewBonds > 0) {
#ifdef DEBUG_FMCS
    std::cerr << "6) Seed::grow numErasedNewBonds " << numErasedNewBonds << std::endl;
#endif
    std::vector<NewBond> dirtyNewBonds;
    dirtyNewBonds.reserve(NewBonds.size());
    dirtyNewBonds.swap(NewBonds);
    for (std::vector<NewBond>::const_iterator nbi = dirtyNewBonds.begin();
         nbi != dirtyNewBonds.end(); nbi++) {
      if (NotSet != nbi->BondIdx) {
        NewBonds.push_back(*nbi);
      }
    }
  }
  // add all other from 2^k-1 possible seeds, where k=newBonds.size()
  // if just one new bond, then seed has already been created
  if (NewBonds.size() > 1) {
    if (sizeof(unsigned long long) * 8 < NewBonds.size()) {
      throw std::runtime_error(
          "Max number of new external bonds of a seed more than 64");
    }
#ifdef DEBUG_FMCS
    std::cerr << "7) Seed::grow NewBonds.size() " << NewBonds.size() << std::endl;
#endif
    BitSet maxCompositionValue;
    Composition2N::compute2N(NewBonds.size(), maxCompositionValue);
    maxCompositionValue -= 1;  // 2^N-1
    Composition2N composition(maxCompositionValue, maxCompositionValue);

#ifdef EXCLUDE_WRONG_COMPOSITION
    std::vector<BitSet> failedCombinations;
    BitSet failedCombinationsMask = 0uLL;
#endif
    while (composition.generateNext()) {
      // exclude already processed single external bond combinations
      if (composition.is2Power()) {
        continue;
      }
      if (0 == numErasedNewBonds &&
          composition.getBitSet() == maxCompositionValue) {
        continue;  // exclude already processed all external bonds combination
      }
      // 2N-1
#ifdef EXCLUDE_WRONG_COMPOSITION
      // OPTIMISATION. reduce amount of generated seeds and match calls
      // 2120 instead of 2208 match calls on small test. 43 wrongComp-s, 83
      // rejected
      if (failedCombinationsMask & composition.getBitSet()) {
        // possibly exists in the list
        bool compositionWrong = false;
        for (std::vector<BitSet>::const_iterator failed =
                 failedCombinations.begin();
             failed != failedCombinations.end(); failed++)
          if (*failed == (*failed & composition.getBitSet())) {
            // combination includes failed combination
            compositionWrong = true;
            break;
          }
        if (compositionWrong) {
#ifdef VERBOSE_STATISTICS_ON
          ++mcs.VerboseStatistics.WrongCompositionRejected;
#endif
          continue;
        }
      }
#endif
#ifdef VERBOSE_STATISTICS_ON
      { ++mcs.VerboseStatistics.Seed; }
#endif
      Seed seed;
      seed.createFromParent(this);
      newAtomsSet.clear();

      for (unsigned i = 0; i < NewBonds.size(); i++) {
        if (composition.isSet(i)) {
          const NewBond* nbi = &NewBonds[i];
          unsigned aIdx =
              nbi->EndAtomIdx;   // existed in this parent seed (ring) or -1
          if (NotSet == aIdx) {  // new atom
            if (newAtomsSet.find(nbi->NewAtomIdx) == newAtomsSet.end()) {
              const Atom* end_atom = nbi->NewAtom;
              aIdx = seed.addAtom(end_atom);
              newAtomsSet.insert(nbi->NewAtomIdx);
            }
          }
          const Bond* src_bond = qmol.getBondWithIdx(nbi->BondIdx);
          seed.addBond(src_bond);
        }
      }
      seed.computeRemainingSize(qmol);
      if (!seed.canGrowBiggerThan(
              mcs.getMaxNumberBonds(),
              mcs.getMaxNumberAtoms())) {  // prune(). // seed too small
#ifdef VERBOSE_STATISTICS_ON
        ++mcs.VerboseStatistics.RemainingSizeRejected;
#endif
      } else {
        seed.MatchResult = MatchResult;
        bool found = mcs.checkIfMatchAndAppend(seed);

        if (!found) {
#ifdef EXCLUDE_WRONG_COMPOSITION  // if seed does not matched it is possible to
                                  // exclude this mismatched combination for
                                  // performance improvement
          failedCombinations.push_back(composition.getBitSet());
          failedCombinationsMask &= composition.getBitSet();
#ifdef VERBOSE_STATISTICS_ON
          ++mcs.VerboseStatistics.WrongCompositionDetected;
#endif
#endif
        }
      }
    }
  }
#ifdef DEBUG_FMCS
  std::cerr << "8) Seed::grow here" << std::endl;
#endif
  GrowingStage = NotSet;  // finished
}

void Seed::computeRemainingSize(const ROMol& qmol) {
  RemainingBonds = RemainingAtoms = 0;

  std::vector<unsigned int> end_atom_stack;
  auto visitedBonds = ExcludedBonds;
  boost::dynamic_bitset<> visitedAtoms(qmol.getNumAtoms());

  std::for_each(MoleculeFragment.AtomsIdx.begin(), MoleculeFragment.AtomsIdx.end(), [&visitedAtoms](unsigned int i) {
    visitedAtoms.set(i);
  });

  // SDF all paths
  // 1. direct neighbours
  for (unsigned seedAtomIdx = LastAddedAtomsBeginIdx;
       seedAtomIdx < getNumAtoms();
       seedAtomIdx++) {  // just now added new border vertices (candidates for
                         // future growing)
    const Atom* atom = MoleculeFragment.Atoms[seedAtomIdx];
    ROMol::OEDGE_ITER beg, end;
    for (boost::tie(beg, end) = qmol.getAtomBonds(atom); beg != end;
         beg++) {  // all bonds from MoleculeFragment.Atoms[srcAtomIdx]
      const Bond& bond = *(qmol[*beg]);
      if (!visitedBonds.test(bond.getIdx())) {
        ++RemainingBonds;
        visitedBonds.set(bond.getIdx());
        unsigned end_atom_idx =
            (MoleculeFragment.AtomsIdx[seedAtomIdx] == bond.getBeginAtomIdx())
                ? bond.getEndAtomIdx()
                : bond.getBeginAtomIdx();
        if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
          ++RemainingAtoms;
          visitedAtoms.set(end_atom_idx);
          end_atom_stack.push_back(end_atom_idx);
        }
      }
    }
  }
  // 2. go deep
  while (!end_atom_stack.empty()) {
    unsigned ai = end_atom_stack.back();
    end_atom_stack.pop_back();
    const Atom* atom = qmol.getAtomWithIdx(ai);
    ROMol::OEDGE_ITER beg, end;
    for (boost::tie(beg, end) = qmol.getAtomBonds(atom); beg != end;
         beg++) {  // all bonds from end_atom
      const Bond& bond = *(qmol[*beg]);
      if (!visitedBonds.test(bond.getIdx())) {
        ++RemainingBonds;
        visitedBonds.set(bond.getIdx());
        unsigned end_atom_idx = (ai == bond.getBeginAtomIdx())
                                    ? bond.getEndAtomIdx()
                                    : bond.getBeginAtomIdx();
        if (!visitedAtoms.test(end_atom_idx)) {  // check RING/CYCLE
          ++RemainingAtoms;
          visitedAtoms.set(end_atom_idx);
          end_atom_stack.push_back(end_atom_idx);
        }
      }
    }
  }
}
}  // namespace FMCS
}  // namespace RDKit
