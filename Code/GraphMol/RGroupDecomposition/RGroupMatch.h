//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_MATCH_DATA
#define RGROUP_MATCH_DATA
#include "RGroupData.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
typedef boost::shared_ptr<RGroupData> RData;
typedef std::map<int, RData> R_DECOMP;

//! RGroupMatch is the decomposition for a single molecule
struct RGroupMatch {
 public:
  size_t core_idx;  // index of the matching core
  size_t numberMissingUserRGroups;
  R_DECOMP rgroups;        // rlabel->RGroupData mapping
  RWMOL_SPTR matchedCore;  // Core with dummy or query atoms and bonds matched

  RGroupMatch(size_t core_index, size_t numberMissingUserRGroups,
              R_DECOMP input_rgroups, RWMOL_SPTR matchedCore)
      : core_idx(core_index),
        numberMissingUserRGroups(numberMissingUserRGroups),
        rgroups(std::move(input_rgroups)),
        matchedCore(std::move(matchedCore)) {}

  std::string toString() const {
    auto rGroupsString = std::accumulate(
        rgroups.cbegin(), rgroups.cend(), std::string(),
        [](std::string s, const std::pair<int, RData> &rgroup) {
          return std::move(s) + "\n\t(" + std::to_string(rgroup.first) + ':' +
                 rgroup.second->toString() + ')';
        });
    std::stringstream ss;
    ss << "Match coreIdx " << core_idx << " missing count "
       << numberMissingUserRGroups << " " << rGroupsString;
    return ss.str();
  }

  void setInputMoleculeForHighlights(const RWMOL_SPTR &inputMol) {
    inputMolForHighlights = inputMol;
    inputMolWasTrimmed = false;
  }

  RWMOL_SPTR getInputMoleculeForHighlights(bool trimHs) {
    if (!inputMolForHighlights || !trimHs || inputMolWasTrimmed) {
      return inputMolForHighlights;
    }
    int numMolAtoms = inputMolForHighlights->getNumAtoms();
    std::vector<int> storedAtomMapNums(numMolAtoms);
    std::vector<std::pair<int, int>> oldBondEnds(
        inputMolForHighlights->getNumBonds(), std::make_pair(-1, -1));
    auto atoms = inputMolForHighlights->atoms();
    std::transform(atoms.begin(), atoms.end(), storedAtomMapNums.begin(),
                   [](auto atom) {
                     auto res = atom->getAtomMapNum();
                     atom->setAtomMapNum(0);
                     return res;
                   });
    for (const auto &pair : rgroups) {
      auto &combinedMol = pair.second->combinedMol;
      std::vector<int> bondIndices;
      if (combinedMol->getPropIfPresent(common_properties::_rgroupTargetBonds,
                                        bondIndices)) {
        std::for_each(bondIndices.begin(), bondIndices.end(),
                      [this, &oldBondEnds](const auto &bondIdx) {
                        const auto bond =
                            inputMolForHighlights->getBondWithIdx(bondIdx);
                        CHECK_INVARIANT(bond, "bond must not be null");
                        const auto beginAtom = bond->getBeginAtom();
                        const auto endAtom = bond->getEndAtom();
                        oldBondEnds[bondIdx].first = beginAtom->getIdx();
                        oldBondEnds[bondIdx].second = endAtom->getIdx();
                        beginAtom->setAtomMapNum(beginAtom->getIdx() + 1);
                        endAtom->setAtomMapNum(endAtom->getIdx() + 1);
                      });
      }
    }
    std::vector<int> oldToNewAtomIndices(numMolAtoms, -1);
    MolOps::RemoveHsParameters rhps;
    rhps.removeMapped = false;
    MolOps::removeHs(*inputMolForHighlights, rhps);
    for (auto atom : inputMolForHighlights->atoms()) {
      auto atomMapNum = atom->getAtomMapNum();
      if (atomMapNum) {
        --atomMapNum;
        oldToNewAtomIndices[atomMapNum] = atom->getIdx();
        atom->setAtomMapNum(storedAtomMapNums.at(atomMapNum));
      }
    }
    for (const auto &pair : rgroups) {
      auto &combinedMol = pair.second->combinedMol;
      std::vector<int> atomIndices;
      if (combinedMol->getPropIfPresent(common_properties::_rgroupTargetAtoms,
                                        atomIndices)) {
        std::transform(
            atomIndices.begin(), atomIndices.end(), atomIndices.begin(),
            [&oldToNewAtomIndices](auto &atomIdx) {
              auto newAtomIdx = oldToNewAtomIndices.at(atomIdx);
              CHECK_INVARIANT(newAtomIdx != -1, "newAtomIdx must be >=0");
              return newAtomIdx;
            });
      }
      combinedMol->setProp(common_properties::_rgroupTargetAtoms, atomIndices);
      std::vector<int> bondIndices;
      if (combinedMol->getPropIfPresent(common_properties::_rgroupTargetBonds,
                                        bondIndices)) {
        std::transform(
            bondIndices.begin(), bondIndices.end(), bondIndices.begin(),
            [this, &oldBondEnds, &oldToNewAtomIndices](auto &bondIdx) {
              const auto &oldPair = oldBondEnds.at(bondIdx);
              CHECK_INVARIANT(oldPair.first != -1 && oldPair.second != -1,
                              "oldPair members must be >=0");
              const auto newBeginAtomIdx =
                  oldToNewAtomIndices.at(oldPair.first);
              const auto newEndAtomIdx = oldToNewAtomIndices.at(oldPair.second);
              CHECK_INVARIANT(newBeginAtomIdx != -1 && newEndAtomIdx != -1,
                              "newBeginAtomIdx and newEndAtomIdx must be >=0");
              const auto bond = inputMolForHighlights->getBondBetweenAtoms(
                  newBeginAtomIdx, newEndAtomIdx);
              CHECK_INVARIANT(bond, "bond must not be null");
              return bond->getIdx();
            });
      }
      combinedMol->setProp(common_properties::_rgroupTargetBonds, bondIndices);
    }
    inputMolWasTrimmed = true;
    return inputMolForHighlights;
  }

 private:
  bool inputMolWasTrimmed = false;
  RWMOL_SPTR inputMolForHighlights;
};

}  // namespace RDKit
#endif
