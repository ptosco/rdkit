//
//  Copyright (c) 2017-2023, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "RGroupData.h"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <regex>

namespace RDKit {

void RGroupData::updateHighlightMaps(ROMOL_SPTR mol, const std::vector<int> &rlabel_attachments) {
  CHECK_INVARIANT(combinedMol, "combinedMol must not be null");
  std::vector<int> atomIndices;
  std::vector<int> bondIndices;
  bool hasAtomIndices = mol->getPropIfPresent(common_properties::_rgroupTargetAtoms, atomIndices);
  bool hasBondIndices = mol->getPropIfPresent(common_properties::_rgroupTargetBonds, bondIndices);
  if (hasAtomIndices) {
    std::transform(rlabel_attachments.begin(), rlabel_attachments.end(), std::inserter(rlabelAtomIndicesMap, rlabelAtomIndicesMap.end()), [&atomIndices](int rlabel) {
      return std::make_pair(rlabel, atomIndices);
    });
  }
  if (hasBondIndices) {
    std::transform(rlabel_attachments.begin(), rlabel_attachments.end(), std::inserter(rlabelBondIndicesMap, rlabelBondIndicesMap.end()), [&bondIndices](int rlabel) {
      return std::make_pair(rlabel, bondIndices);
    });
  }
}

void RGroupData::setRGroupHighlightsAsJSON(const std::map<int, int> &rlabelMapping) {
  combinedMol->setProp(common_properties::_rgroupTargetAtoms, getHighlightsAsJSON(rlabelMapping, rlabelAtomIndicesMap));
  combinedMol->setProp(common_properties::_rgroupTargetBonds, getHighlightsAsJSON(rlabelMapping, rlabelBondIndicesMap));
}

std::string RGroupData::getHighlightsAsJSON(const std::map<int, int> &rlabelMapping, const std::map<int, std::vector<int>> &rlabelIndicesMap) {
  static const std::string QUOTE{"\""};
  static const std::string QUOTE_COLON_OPEN_SQUARE_BRACKET{"\":["};
  static const std::string CLOSED_SQUARE_BRACKET{"]"};
  static const std::string COMMA{","};
  static const std::string OPEN_CURLY_BRACKET{"{"};
  static const std::string CLOSED_CURLY_BRACKET{"}"};
  std::string res = OPEN_CURLY_BRACKET;
  auto rlabelIndicesMapSize = rlabelIndicesMap.size();
  for (const auto& [rlabel, indices] : rlabelIndicesMap) {
    auto rlabelTmpToFinalPair = rlabelMapping.find(rlabel);
    if (rlabelTmpToFinalPair == rlabelMapping.end()) {
      continue;
    }
    res += QUOTE + getRGroupLabel(rlabelTmpToFinalPair->second) + QUOTE_COLON_OPEN_SQUARE_BRACKET;
    for (auto indicesIt = indices.begin(); indicesIt != indices.end(); ++indicesIt) {
      res += std::to_string(*indicesIt) + (indicesIt == indices.end() - 1 ? CLOSED_SQUARE_BRACKET : COMMA);
    }
    res += (--rlabelIndicesMapSize ? COMMA : CLOSED_CURLY_BRACKET);
  }
  return res;
}

std::string RGroupData::getRGroupLabel(int rlabel) {
  static const std::string RPREFIX = "R";
  return RPREFIX + std::to_string(rlabel);
}

const std::string &RGroupData::getCoreLabel() {
  static const std::string CORE = "Core";
  return CORE;
}

void RGroupData::add(const ROMOL_SPTR &newMol,
                     const std::vector<int> &rlabel_attachments) {
  // some fragments can be added multiple times if they are cyclic
  if (std::any_of(mols.begin(), mols.end(), [&newMol](const auto &mol) {
    return newMol == mol;
  })) {
    return;
  }

  if (!mols.empty()) {
    // don't add extraneous hydrogens
    if (isMolHydrogen(*newMol)) {
      return;
    }
    if (is_hydrogen) {
      // if we are adding a heavy attachment to hydrogens, discard the
      // hydrogen and start over
      combinedMol = nullptr;
      smilesVect.clear();
      attachments.clear();
      mols.clear();
    }
  }

  labelled = false;
  std::copy(rlabel_attachments.begin(), rlabel_attachments.end(),
            std::inserter(attachments, attachments.end()));

  mols.push_back(newMol);
  static const std::regex remove_isotopes_regex("\\[\\d*\\*\\]");
  // remove the isotope labels from the SMILES string to avoid
  // that identical R-group are perceived as different when
  // MCS alignment is not used (NoAlign flag)
  smilesVect.push_back(std::regex_replace(MolToSmiles(*newMol, true),
                                          remove_isotopes_regex, "*"));
  if (!combinedMol) {
    combinedMol = RWMOL_SPTR(new RWMol(*newMol));
  } else {
    combinedMol.reset(static_cast<RWMol *>(combineMols(*combinedMol, *newMol)));
    single_fragment = false;
  }
  updateHighlightMaps(newMol, rlabel_attachments);
  smiles = getSmiles();
  combinedMol->setProp(common_properties::internalRgroupSmiles, smiles);
  computeIsHydrogen();
  is_linker = single_fragment && attachments.size() > 1;
}

std::map<int, int> RGroupData::getNumBondsToRlabels() const {
  std::map<int, int> rlabelsUsedCount;

  for (const auto atom : combinedMol->atoms()) {
    int rlabel;
    if (atom->getPropIfPresent<int>(RLABEL, rlabel)) {
      ++rlabelsUsedCount[rlabel];
    }
  }
  return rlabelsUsedCount;
}

std::string RGroupData::toString() const {
  auto attachmentString = std::accumulate(
      attachments.cbegin(), attachments.cend(), std::string(),
      [](std::string s, int a) {
        return s.empty() ? std::to_string(a)
                         : std::move(s) + ',' + std::to_string(a);
      });
  std::stringstream ss;
  ss << "RG " << attachmentString << " " << getSmiles();
  return ss.str();
}

void RGroupData::computeIsHydrogen() {  // is the rgroup all Hs
  is_hydrogen = std::all_of(mols.begin(), mols.end(), [this](const auto &mol) {
    return isMolHydrogen(*mol);
  });
}

bool RGroupData::isMolHydrogen(const ROMol &mol) const {
  auto atoms = mol.atoms();
  return std::all_of(atoms.begin(), atoms.end(), [](const auto &atom) {
    return (atom->getAtomicNum() == 1 || (atom->getAtomicNum() == 0 && atom->hasProp(SIDECHAIN_RLABELS)));
  });
}

//! compute the canonical smiles for the attachments (bug: removes dupes since
//! we are using a set...)
std::string RGroupData::getSmiles() const {
  std::string s;
  for (const auto &it : smilesVect) {
    if (s.length()) {
      s += ".";
    }
    s += it;
  }
  return s;
}

}  // namespace RDKit
