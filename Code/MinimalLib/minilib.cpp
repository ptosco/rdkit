//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string>
#include <iostream>
#include "minilib.h"

#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolDraw2D/MolDraw2D.h>
#include <GraphMol/MolDraw2D/MolDraw2DSVG.h>
#include <GraphMol/MolDraw2D/MolDraw2DUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <DataStructs/BitOps.h>

#include <INCHI-API/inchi.h>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include <boost/dynamic_bitset.hpp>

namespace rj = rapidjson;

using namespace RDKit;

std::string process_details(const std::string &details, unsigned int &width,
                            unsigned int &height, int &offsetx, int &offsety,
                            std::string &legend, std::vector<int> &atomIds,
                            std::vector<int> &bondIds) {
  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) return "Invalid JSON";

  if (doc.HasMember("atoms")) {
    if (!doc["atoms"].IsArray()) {
      return "JSON doesn't contain 'atoms' field, or it is not an array";
    }
    for (const auto &molval : doc["atoms"].GetArray()) {
      if (!molval.IsInt()) return ("Atom IDs should be integers");
      atomIds.push_back(molval.GetInt());
    }
  }
  if (doc.HasMember("bonds")) {
    if (!doc["bonds"].IsArray()) {
      return "JSON contain 'bonds' field, but it is not an array";
    }
    for (const auto &molval : doc["bonds"].GetArray()) {
      if (!molval.IsInt()) return ("Bond IDs should be integers");
      bondIds.push_back(molval.GetInt());
    }
  }

  if (doc.HasMember("width")) {
    if (!doc["width"].IsUint()) {
      return "JSON contains 'width' field, but it is not an unsigned int";
    }
    width = doc["width"].GetUint();
  }

  if (doc.HasMember("height")) {
    if (!doc["height"].IsUint()) {
      return "JSON contains 'height' field, but it is not an unsigned int";
    }
    height = doc["height"].GetUint();
  }

  if (doc.HasMember("offsetx")) {
    if (!doc["offsetx"].IsInt()) {
      return "JSON contains 'offsetx' field, but it is not an int";
    }
    offsetx = doc["offsetx"].GetInt();
  }

  if (doc.HasMember("offsety")) {
    if (!doc["offsety"].IsInt()) {
      return "JSON contains 'offsety' field, but it is not an int";
    }
    offsety = doc["offsety"].GetInt();
  }

  if (doc.HasMember("legend")) {
    if (!doc["legend"].IsString()) {
      return "JSON contains 'legend' field, but it is not a string";
    }
    legend = doc["legend"].GetString();
  }

  return "";
}

namespace {
RWMol *mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else {
    SmilesParserParams ps;
    ps.sanitize = false;
    res = SmilesToMol(input, ps);
  }
  if (res) {
    try {
      MolOps::sanitizeMol(*res);
      MolOps::assignStereochemistry(*res, true, true, true);
    } catch (...) {
      delete res;
      res = nullptr;
    }
  }
  return res;
}

RWMol *qmol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else {
    res = SmartsToMol(input);
  }
  return res;
}

std::string svg_(const ROMol &m, unsigned int w, unsigned int h,
                 const std::string &details = "") {
  std::vector<int> atomIds;
  std::vector<int> bondIds;
  std::string legend = "";
  int offsetx = 0, offsety = 0;
  if (!details.empty()) {
    auto problems = process_details(details, w, h, offsetx, offsety, legend,
                                    atomIds, bondIds);
    if (!problems.empty()) {
      return problems;
    }
  }

  MolDraw2DSVG drawer(w, h);
  if (!details.empty()) {
    MolDraw2DUtils::updateDrawerParamsFromJSON(drawer, details);
  }
  drawer.setOffset(offsetx, offsety);

  MolDraw2DUtils::prepareAndDrawMolecule(drawer, m, legend, &atomIds, &bondIds);
  drawer.finishDrawing();

  return drawer.getDrawingText();
}
}  // namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return MolToSmiles(*d_mol);
}
std::string JSMol::get_cxsmiles() const {
  if (!d_mol) return "";
  return MolToCXSmiles(*d_mol);
}
std::string JSMol::get_svg(unsigned int w, unsigned int h) const {
  if (!d_mol) return "";
  return svg_(*d_mol, w, h);
}
std::string JSMol::get_svg_with_highlights(const std::string &details) const {
  if (!d_mol) return "";

  unsigned int w = d_defaultWidth;
  unsigned int h = d_defaultHeight;
  return svg_(*d_mol, w, h, details);
}

std::string JSMol::get_inchi() const {
  if (!d_mol) return "";
  ExtraInchiReturnValues rv;
  return MolToInchi(*d_mol, rv);
}
std::string JSMol::get_molblock() const {
  if (!d_mol) return "";
  return MolToMolBlock(*d_mol);
}
std::string JSMol::get_v3Kmolblock() const {
  if (!d_mol) return "";
  return MolToV3KMolBlock(*d_mol);
}

namespace {
void get_sss_json(const ROMol *d_mol, const ROMol *q_mol,
                  const MatchVectType &match, rj::Value &obj,
                  rj::Document &doc) {
  rj::Value rjAtoms(rj::kArrayType);
  for (const auto &pr : match) {
    rjAtoms.PushBack(pr.second, doc.GetAllocator());
  }
  obj.AddMember("atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto qbond : q_mol->bonds()) {
    unsigned int idx1 = match[qbond->getBeginAtomIdx()].second;
    unsigned int idx2 = match[qbond->getEndAtomIdx()].second;
    const auto bond = d_mol->getBondBetweenAtoms(idx1, idx2);
    if (bond != nullptr) {
      rjBonds.PushBack(bond->getIdx(), doc.GetAllocator());
    }
  }
  obj.AddMember("bonds", rjBonds, doc.GetAllocator());
}
}  // namespace

std::string JSMol::get_substruct_match(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) return res;

  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    rj::Document doc;
    doc.SetObject();
    get_sss_json(d_mol.get(), q.d_mol.get(), match, doc, doc);
    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return res;
}

std::string JSMol::get_substruct_matches(const JSMol &q) const {
  std::string res = "{}";
  if (!d_mol || !q.d_mol) return res;

  auto matches = SubstructMatch(*d_mol, (*q.d_mol));
  if (!matches.empty()) {
    rj::Document doc;
    doc.SetArray();

    for (const auto &match : matches) {
      rj::Value rjMatch(rj::kObjectType);
      get_sss_json(d_mol.get(), q.d_mol.get(), match, rjMatch, doc);
      doc.PushBack(rjMatch, doc.GetAllocator());
    }

    rj::StringBuffer buffer;
    rj::Writer<rj::StringBuffer> writer(buffer);
    doc.Accept(writer);
    res = buffer.GetString();
  }

  return res;
}

std::string JSMol::get_descriptors() const {
  if (!d_mol) return "{}";
  rj::Document doc;
  doc.SetObject();

  Descriptors::Properties props;
  std::vector<std::string> dns = props.getPropertyNames();
  std::vector<double> dvs = props.computeProperties(*d_mol);
  for (size_t i = 0; i < dns.size(); ++i) {
    rj::Value v(dvs[i]);
    const auto srt = rj::StringRef(dns[i].c_str());
    doc.AddMember(srt, v, doc.GetAllocator());
  }

  if (std::find(dns.begin(), dns.end(), std::string("amw")) == dns.end()) {
    rj::Value v(Descriptors::calcAMW(*d_mol));
    doc.AddMember("amw", v, doc.GetAllocator());
  }

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  writer.SetMaxDecimalPlaces(5);
  doc.Accept(writer);
  return buffer.GetString();
}

std::string JSMol::get_morgan_fp(unsigned int radius,
                                 unsigned int fplen) const {
  if (!d_mol) return "";
  auto fp = MorganFingerprints::getFingerprintAsBitVect(*d_mol, radius, fplen);
  std::string res = BitVectToText(*fp);
  delete fp;
  return res;
}

std::string JSMol::get_stereo_tags() const {
  if (!d_mol) return "{}";
  rj::Document doc;
  doc.SetObject();

  bool cleanIt = true;
  bool force = true;
  bool flagPossibleStereocenters = true;
  MolOps::assignStereochemistry(*d_mol, cleanIt, force,
                                flagPossibleStereocenters);
  CIPLabeler::assignCIPLabels(*d_mol);

  rj::Value rjAtoms(rj::kArrayType);
  for (const auto atom : d_mol->atoms()) {
    std::string cip;
    if (!atom->getPropIfPresent(common_properties::_CIPCode, cip)) {
      if (atom->hasProp(common_properties::_ChiralityPossible)) {
        cip = "?";
      }
    }
    if (!cip.empty()) {
      cip = "(" + cip + ")";
      rj::Value entry(rj::kArrayType);
      entry.PushBack(atom->getIdx(), doc.GetAllocator());
      rj::Value v;
      v.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
      entry.PushBack(v, doc.GetAllocator());
      rjAtoms.PushBack(entry, doc.GetAllocator());
    }
  }
  doc.AddMember("CIP_atoms", rjAtoms, doc.GetAllocator());

  rj::Value rjBonds(rj::kArrayType);
  for (const auto bond : d_mol->bonds()) {
    std::string cip;
    if (bond->getPropIfPresent(common_properties::_CIPCode, cip)) {
      cip = "(" + cip + ")";
      rj::Value entry(rj::kArrayType);
      entry.PushBack(bond->getBeginAtomIdx(), doc.GetAllocator());
      entry.PushBack(bond->getEndAtomIdx(), doc.GetAllocator());
      rj::Value v;
      v.SetString(cip.c_str(), cip.size(), doc.GetAllocator());
      entry.PushBack(v, doc.GetAllocator());
      rjBonds.PushBack(entry, doc.GetAllocator());
    }
  }

  doc.AddMember("CIP_bonds", rjBonds, doc.GetAllocator());

  rj::StringBuffer buffer;
  rj::Writer<rj::StringBuffer> writer(buffer);
  doc.Accept(writer);
  return buffer.GetString();
}

std::string JSMol::get_aromatic_form() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::setAromaticity(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = false;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::get_kekule_form() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::Kekulize(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::get_new_coords(bool useCoordGen) const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif
  RDDepict::compute2DCoords(molCopy);
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif

  return MolToMolBlock(molCopy);
}

std::string JSMol::remove_hs() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::removeAllHs(molCopy);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::add_hs() const {
  if (!d_mol) return "";

  RWMol molCopy(*d_mol);
  MolOps::addHs(molCopy);

  // RDDepict::generateDepictionMatching2DStructure(molCopy, *d_mol);

  bool includeStereo = true;
  int confId = -1;
  bool kekulize = true;
  return MolToMolBlock(molCopy, includeStereo, confId, kekulize);
}

std::string JSMol::condense_abbreviations(double maxCoverage, bool useLinkers) {
  if (!d_mol) return "";
  if (!useLinkers) {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultAbbreviations(), maxCoverage);
  } else {
    Abbreviations::condenseMolAbbreviations(
        *d_mol, Abbreviations::Utils::getDefaultLinkers(), maxCoverage);
  }
  return "";
}

std::string JSMol::condense_abbreviations_from_defs(
    const std::string &definitions, double maxCoverage, bool areLinkers) {
  static std::string lastDefs = "";
  static std::vector<Abbreviations::AbbreviationDefinition> abbrevs;
  if (definitions != lastDefs) {
    // yes, we are making the assumption that the "areLinkers" argument remains
    // the same if the definitions are the same
    bool removeExtraDummies = areLinkers;
    bool allowConnectionToDummies = areLinkers;
    lastDefs = definitions;
    try {
      abbrevs = Abbreviations::Utils::parseAbbreviations(
          definitions, removeExtraDummies, allowConnectionToDummies);
    } catch (...) {
      return "cannot parse abbreviations";
    }
  }
  Abbreviations::condenseMolAbbreviations(*d_mol, abbrevs, maxCoverage);
}

class MatchScorer {
 public:
  MatchScorer(const ROMol &mol, const ROMol &query) : d_mol(mol), d_query(query), d_sumIndices(0) {
    for (auto i = 0; i < d_mol.getNumAtoms(); ++i) {
      d_sumIndices += static_cast<double>(i);
    }
  }
  double score(const MatchVectType &match) {
    auto it = d_cache.find(&match);
    if (it == d_cache.end()) {
      it = computeAndCacheScore(match);
    }
    return it->second;
  }
  bool doesQueryRGroupMatchHydrogen(const std::pair<int, int> &pair) {
    const auto queryAtom = d_query.getAtomWithIdx(pair.first);
    const auto molAtom = d_mol.getAtomWithIdx(pair.second);
    return (queryAtom->hasQuery() && queryAtom->getDegree() == 1 && molAtom->getAtomicNum() == 1);
  }
 private:
  std::unordered_map<const MatchVectType *, double>::iterator computeAndCacheScore(const MatchVectType &match) {
    double penalty = 0.0;
    double i = 0.0;
    for (const auto &pair : match) {
      i += static_cast<double>(pair.second);
      if (doesQueryRGroupMatchHydrogen(pair)) {
        penalty += 1.0;
      }
    }
    penalty += i / d_sumIndices;
    return d_cache.insert(std::make_pair(&match, penalty)).first;
  }
  std::unordered_map<const MatchVectType *, double> d_cache;
  const ROMol &d_mol;
  const ROMol &d_query;
  double d_sumIndices;
};

std::string JSMol::generate_aligned_coords(const JSMol &templateMol, bool useCoordGen, bool allowRGroups){
  if (!d_mol || !templateMol.d_mol || !templateMol.d_mol->getNumConformers()) return "";

#ifdef RDK_BUILD_COORDGEN_SUPPORT
  bool oprefer = RDDepict::preferCoordGen;
  RDDepict::preferCoordGen = useCoordGen;
#endif 
  RDKit::ROMol *refPattern = nullptr;
  std::unique_ptr<RDKit::RWMol> templateCopy;
  RDKit::RWMol *templateRWMol = templateMol.d_mol.get();
  if (allowRGroups) {
    std::unique_ptr<RDKit::ROMol> molHs(RDKit::MolOps::addHs(*static_cast<ROMol *>(d_mol.get())));
    if (molHs) {
      MatchScorer matchScorer(*molHs, *templateRWMol);
      auto matches = SubstructMatch(*molHs, *templateRWMol);
      auto bestMatch = *std::min_element(matches.begin(), matches.end(),
        [&matchScorer](const MatchVectType &aMatch, const MatchVectType &bMatch) {
          return (matchScorer.score(aMatch) < matchScorer.score(bMatch));
        });
      if (matchScorer.score(bestMatch) > 1.0) {
        templateCopy.reset(new RDKit::RWMol(*templateRWMol));
        templateRWMol = templateCopy.get();
        for (const auto &pair : bestMatch) {
          if (matchScorer.doesQueryRGroupMatchHydrogen(pair)) {
            Atom *dummyAtom = templateRWMol->getAtomWithIdx(pair.first);
            dummyAtom->setAtomicNum(1);
            dummyAtom->setIsotope(0);
          }
        }
        RDKit::MolOps::removeAllHs(*templateRWMol);
      }
    }
  }
  bool acceptFailure = true;
  int confId = -1;
  RDDepict::generateDepictionMatching2DStructure(*d_mol, *templateRWMol, confId,
     refPattern, acceptFailure);
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = oprefer;
#endif
  return "";
}


std::string get_inchikey_for_inchi(const std::string &input) {
  return InchiToInchiKey(input);
}

JSMol *get_mol(const std::string &input) {
  RWMol *mol = mol_from_input(input);
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  RWMol *mol = qmol_from_input(input);
  return new JSMol(mol);
}

std::string version() { return std::string(rdkitVersion); }

void prefer_coordgen(bool useCoordGen) {
#ifdef RDK_BUILD_COORDGEN_SUPPORT
  RDDepict::preferCoordGen = useCoordGen;
#endif
}
