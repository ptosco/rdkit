// $Id$
//
//  Copyright (C) 2013-2016 Paolo Tosco
//
//  Copyright (C) 2004-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <iostream>
#include <cmath>
#include <cctype>

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/Contribs.h>
#include "AtomTyper.h"
#include "Builder.h"
#ifdef RDK_BUILD_WITH_OPENMM
#include <OpenMM.h>
#include <OpenMMMMFF.h>
#include <openmm/Units.h>
#endif

//#define PRINT_MMFF_FORCES 1

namespace RDKit {
namespace MMFF {
using namespace ForceFields::MMFF;

namespace Tools {
// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
#ifdef RDK_BUILD_WITH_OPENMM
class OpenMMPlugins {
  public:
    const std::vector<std::string>& load(
      bool force = false, const std::string &dir = std::string());
    const std::vector<std::string>& loadedPlugins() {
      return d_loadedPlugins;
    }
    const std::vector<std::string>& failedPlugins() {
      return d_failedPlugins;
    }
    static OpenMMPlugins *instance();
  private:
    //! to force this to be a singleton, the constructor must be private
    OpenMMPlugins() :
      d_alreadyLoaded(false) {};
    static class OpenMMPlugins *ds_instance;  //!< the singleton
    bool d_alreadyLoaded;
    std::vector<std::string> d_loadedPlugins;
    std::vector<std::string> d_failedPlugins;
};

OpenMMPlugins *OpenMMPlugins::ds_instance = NULL;

OpenMMPlugins *OpenMMPlugins::instance()
{
  if (!ds_instance)
    ds_instance = new OpenMMPlugins();
  return ds_instance;
}

const std::vector<std::string> &OpenMMPlugins::load(
  bool force, const std::string &dir) {
  if (!(d_alreadyLoaded || force)) {
    std::string d(dir.length() ? dir
      : OpenMM::Platform::getDefaultPluginsDirectory());
    d_loadedPlugins = OpenMM::Platform::loadPluginsFromDirectory(d);
    d_failedPlugins = OpenMM::Platform::getPluginLoadFailures();
    d_alreadyLoaded = true;
  }
  return d_loadedPlugins;
}
#endif

void checkFFPreconditions(ForceFields::ForceField *field,
                          MMFFMolProperties *mmffMolProperties, int ffOpts) {
#ifdef RDK_BUILD_WITH_OPENMM
  RDUNUSED_PARAM(ffOpts);
#else
  PRECONDITION(!(ffOpts & ForceFields::USE_OPENMM),
               "OpenMM support was not built in the RDKit");
#endif
  PRECONDITION(field, "bad ForceField");
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");
}

OpenMMForceField *getOpenMMForceField(ForceFields::ForceField *field, int ffOpts) {
  OpenMMForceField *res = NULL;
#ifdef RDK_BUILD_WITH_OPENMM
  if (ffOpts & ForceFields::USE_OPENMM) {
    res = static_cast<OpenMMForceField *>(field);
  }
#endif
  return res;
}

void addBonds(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
              ForceFields::ForceField *field) {
  addBonds(mol, mmffMolProperties, field, 0);
}

void addBonds(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
              ForceFields::ForceField *field, int ffOpts) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ffOpts);
  double totalBondStretchEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "B O N D   S T R E T C H I N G\n\n"
                 "------ATOMS------   ATOM TYPES   FF     BOND     IDEAL       "
                 "                 FORCE\n"
                 "  I        J          I    J   CLASS   LENGTH   LENGTH    "
                 "DIFF.    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "------------------------" << std::endl;
    }
  }
  for (ROMol::ConstBondIterator bi = mol.beginBonds(); bi != mol.endBonds();
       ++bi) {
    unsigned int idx1 = (*bi)->getBeginAtomIdx();
    unsigned int idx2 = (*bi)->getEndAtomIdx();
    unsigned int bondType;
    MMFFBond mmffBondParams;
    if (mmffMolProperties->getMMFFBondStretchParams(mol, idx1, idx2, bondType,
                                                    mmffBondParams)) {
      if (!fieldOMM) {
        BondStretchContrib *contrib = new BondStretchContrib(field, idx1, idx2, &mmffBondParams);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
#ifdef RDK_BUILD_WITH_OPENMM
      else {
        fieldOMM->addBondStretchContrib(idx1, idx2, &mmffBondParams);
      }
#endif
      if (mmffMolProperties->getMMFFVerbosity()) {
        unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx1);
        unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx2);
        const Atom *iAtom = (*bi)->getBeginAtom();
        const Atom *jAtom = (*bi)->getEndAtom();
        const double dist = field->distance(idx1, idx2);
        const double bondStretchEnergy = MMFF::Utils::calcBondStretchEnergy(
            mmffBondParams.r0, mmffBondParams.kb, dist);
        if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
          oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                  << std::setw(5) << idx1 + 1 << std::setw(2)
                  << jAtom->getSymbol() << " #" << std::setw(5) << idx2 + 1
                  << std::right << std::setw(5) << iAtomType << std::setw(5)
                  << jAtomType << std::setw(6) << bondType << "  " << std::fixed
                  << std::setprecision(3) << std::setw(9) << dist
                  << std::setw(9) << mmffBondParams.r0 << std::setw(9)
                  << dist - mmffBondParams.r0 << std::setw(10)
                  << bondStretchEnergy << std::setw(10) << mmffBondParams.kb
                  << std::endl;
        }
        totalBondStretchEnergy += bondStretchEnergy;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL BOND STRETCH ENERGY      =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalBondStretchEnergy
            << std::endl;
  }
}

unsigned int twoBitCellPos(unsigned int nAtoms, int i, int j) {
  if (j < i) std::swap(i, j);

  return i * (nAtoms - 1) + i * (1 - i) / 2 + j;
}

void setTwoBitCell(boost::shared_array<boost::uint8_t> &res, unsigned int pos,
                   boost::uint8_t value) {
  unsigned int twoBitPos = pos / 4;
  unsigned int shift = 2 * (pos % 4);
  boost::uint8_t twoBitMask = 3 << shift;
  res[twoBitPos] = ((res[twoBitPos] & (~twoBitMask)) | (value << shift));
}

boost::uint8_t getTwoBitCell(boost::shared_array<boost::uint8_t> &res,
                             unsigned int pos) {
  unsigned int twoBitPos = pos / 4;
  unsigned int shift = 2 * (pos % 4);
  boost::uint8_t twoBitMask = 3 << shift;

  return ((res[twoBitPos] & twoBitMask) >> shift);
}

// ------------------------------------------------------------------------
//
// the two-bit matrix returned by this contains:
//   0: if atoms i and j are directly connected
//   1: if atoms i and j are connected via an atom
//   2: if atoms i and j are in a 1,4 relationship
//   3: otherwise
//
//  NOTE: the caller is responsible for calling delete []
//  on the result
//
// ------------------------------------------------------------------------
boost::shared_array<boost::uint8_t> buildNeighborMatrix(const ROMol &mol) {
  const boost::uint8_t RELATION_1_X_INIT = RELATION_1_X
    | (RELATION_1_X << 2) | (RELATION_1_X << 4) | (RELATION_1_X << 6);
  unsigned int nAtoms = mol.getNumAtoms();
  unsigned nTwoBitCells = (nAtoms * (nAtoms + 1) - 1) / 8 + 1;
  boost::shared_array<boost::uint8_t> res(new boost::uint8_t[nTwoBitCells]);
  std::memset(res.get(), RELATION_1_X_INIT, nTwoBitCells);
  for (ROMol::ConstBondIterator bondi = mol.beginBonds(); bondi != mol.endBonds(); ++bondi) {
    setTwoBitCell(res, twoBitCellPos(nAtoms, (*bondi)->getBeginAtomIdx(),
                  (*bondi)->getEndAtomIdx()), RELATION_1_2);
    unsigned int bondiBeginAtomIdx = (*bondi)->getBeginAtomIdx();
    unsigned int bondiEndAtomIdx = (*bondi)->getEndAtomIdx();
    for (ROMol::ConstBondIterator bondj = bondi; ++bondj != mol.endBonds();) {
      int idx1 = -1;
      int idx3 = -1;
      unsigned int bondjBeginAtomIdx = (*bondj)->getBeginAtomIdx();
      unsigned int bondjEndAtomIdx = (*bondj)->getEndAtomIdx();
      if (bondiBeginAtomIdx == bondjBeginAtomIdx) {
        idx1 = bondiEndAtomIdx;
        idx3 = bondjEndAtomIdx;
      } else if (bondiBeginAtomIdx == bondjEndAtomIdx) {
        idx1 = bondiEndAtomIdx;
        idx3 = bondjBeginAtomIdx;
      } else if (bondiEndAtomIdx == bondjBeginAtomIdx) {
        idx1 = bondiBeginAtomIdx;
        idx3 = bondjEndAtomIdx;
      } else if (bondiEndAtomIdx == bondjEndAtomIdx) {
        idx1 = bondiBeginAtomIdx;
        idx3 = bondjBeginAtomIdx;
      } else {
        // check if atoms i and j are in a 1,4-relationship
        if ((mol.getBondBetweenAtoms(bondiBeginAtomIdx,
                                     bondjBeginAtomIdx)) &&
            (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiEndAtomIdx,
                                    bondjEndAtomIdx)) == RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiEndAtomIdx, bondjEndAtomIdx),
              RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiBeginAtomIdx,
                                            bondjEndAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiEndAtomIdx,
                                           bondjBeginAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiEndAtomIdx, bondjBeginAtomIdx),
              RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiEndAtomIdx,
                                            bondjBeginAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiBeginAtomIdx,
                                           bondjEndAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiBeginAtomIdx, bondjEndAtomIdx),
              RELATION_1_4);
        } else if ((mol.getBondBetweenAtoms(bondiEndAtomIdx,
                                            bondjEndAtomIdx)) &&
                   (getTwoBitCell(res, twoBitCellPos(nAtoms, bondiBeginAtomIdx,
                                           bondjBeginAtomIdx)) ==
                    RELATION_1_X)) {
          setTwoBitCell(
              res, twoBitCellPos(nAtoms, bondiBeginAtomIdx, bondjBeginAtomIdx),
              RELATION_1_4);
        }
      }
      if (idx1 > -1) {
        setTwoBitCell(res, twoBitCellPos(nAtoms, idx1, idx3), RELATION_1_3);
      }
    }
  }
  return res;
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addAngles(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
              ForceFields::ForceField *field) {
  addAngles(mol, mmffMolProperties, field, 0);
}

void addAngles(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
               ForceFields::ForceField *field, int ffOpts) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ffOpts);
  unsigned int idx[3];
  MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;

  unsigned int nAtoms = mol.getNumAtoms();
  double totalAngleBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "A N G L E   B E N D I N G\n\n"
                 "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
                 "   IDEAL                          FORCE\n"
                 "  I        J        K          I    J    K   CLASS    ANGLE  "
                 "   ANGLE      DIFF.    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------" << std::endl;
    }
    points = field->positions();
  }
  for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
    if (jAtom->getDegree() == 1) {
      continue;
    }
    unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
    const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
    boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
    for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
      const Atom *iAtom = mol[*nbr1Idx].get();
      idx[0] = iAtom->getIdx();
      boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
      for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
        if (nbr2Idx < (nbr1Idx + 1)) {
          continue;
        }
        const Atom *kAtom = mol[*nbr2Idx].get();
        idx[2] = kAtom->getIdx();
        unsigned int angleType;
        MMFFAngle mmffAngleParams;
        if (mmffMolProperties->getMMFFAngleBendParams(
                mol, idx[0], idx[1], idx[2], angleType, mmffAngleParams)) {
          if (!fieldOMM) {
            AngleBendContrib *contrib = new AngleBendContrib(
                field, idx[0], idx[1], idx[2],
                &mmffAngleParams, mmffPropParamsCentralAtom);
            field->contribs().push_back(ForceFields::ContribPtr(contrib));
          }
#ifdef RDK_BUILD_WITH_OPENMM
          else {
            fieldOMM->addAngleBendContrib(idx[0], idx[1], idx[2],
                &mmffAngleParams, mmffPropParamsCentralAtom);
          }
#endif
          if (mmffMolProperties->getMMFFVerbosity()) {
            unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
            unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
            const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                                     (*(points[idx[0]]))[1],
                                     (*(points[idx[0]]))[2]);
            const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                                     (*(points[idx[1]]))[1],
                                     (*(points[idx[1]]))[2]);
            const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                                     (*(points[idx[2]]))[1],
                                     (*(points[idx[2]]))[2]);
            const double cosTheta =
                MMFF::Utils::calcCosTheta(p1, p2, p3, field->distance(idx[0], idx[1]),
                                          field->distance(idx[1], idx[2]));
            const double theta = RAD2DEG * acos(cosTheta);
            const double angleBendEnergy = MMFF::Utils::calcAngleBendEnergy(
                mmffAngleParams.theta0, mmffAngleParams.ka,
                mmffPropParamsCentralAtom->linh, cosTheta);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                      << std::setw(5) << idx[0] + 1 << std::setw(2)
                      << jAtom->getSymbol() << " #" << std::setw(5)
                      << idx[1] + 1 << std::setw(2) << kAtom->getSymbol()
                      << " #" << std::setw(5) << idx[2] + 1 << std::right
                      << std::setw(5) << iAtomType << std::setw(5) << jAtomType
                      << std::setw(5) << kAtomType << std::setw(6) << angleType
                      << "  " << std::fixed << std::setprecision(3)
                      << std::setw(10) << theta << std::setw(10)
                      << mmffAngleParams.theta0 << std::setw(10)
                      << theta - mmffAngleParams.theta0 << std::setw(10)
                      << angleBendEnergy << std::setw(10) << mmffAngleParams.ka
                      << std::endl;
            }
            totalAngleBendEnergy += angleBendEnergy;
          }
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL ANGLE BEND ENERGY        =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalAngleBendEnergy
            << std::endl;
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addStretchBend(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                    ForceFields::ForceField *field) {
  addStretchBend(mol, mmffMolProperties, field, 0);
}

void addStretchBend(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                    ForceFields::ForceField *field, int ffOpts) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ffOpts);
  unsigned int idx[3];
  MMFFPropCollection *mmffProp = MMFFPropCollection::getMMFFProp();
  std::pair<bool, const MMFFStbn *> mmffStbnParams;
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;

  unsigned int nAtoms = mol.getNumAtoms();
  double totalStretchBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "S T R E T C H   B E N D I N G\n\n"
                 "-----------ATOMS-----------    ATOM TYPES      FF    VALENCE "
                 "    DELTA     DELTA     DELTA               F CON\n"
                 "  I        J        K          I    J    K   CLASS    ANGLE  "
                 "    ANGLE     R(I,J)    R(J,K)   ENERGY    I-J (J-K)\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------------------"
              << std::endl;
    }
    points = field->positions();
  }
  for (idx[1] = 0; idx[1] < nAtoms; ++idx[1]) {
    const Atom *jAtom = mol.getAtomWithIdx(idx[1]);
    if (jAtom->getDegree() == 1) {
      continue;
    }
    unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
    const MMFFProp *mmffPropParamsCentralAtom = (*mmffProp)(jAtomType);
    if (mmffPropParamsCentralAtom->linh) {
      continue;
    }
    boost::tie(nbr1Idx, end1Nbrs) = mol.getAtomNeighbors(jAtom);
    unsigned int i = 0;
    for (; nbr1Idx != end1Nbrs; ++nbr1Idx) {
      const Atom *iAtom = mol[*nbr1Idx].get();
      boost::tie(nbr2Idx, end2Nbrs) = mol.getAtomNeighbors(jAtom);
      unsigned int j = 0;
      for (; nbr2Idx != end2Nbrs; ++nbr2Idx) {
        const Atom *kAtom = mol[*nbr2Idx].get();
        if (j < (i + 1)) {
          ++j;
          continue;
        }
        idx[0] = iAtom->getIdx();
        idx[2] = kAtom->getIdx();
        unsigned int stretchBendType;
        MMFFStbn mmffStbnParams;
        MMFFBond mmffBondParams[2];
        MMFFAngle mmffAngleParams;
        if (mmffMolProperties->getMMFFStretchBendParams(
                mol, idx[0], idx[1], idx[2], stretchBendType, mmffStbnParams,
                mmffBondParams, mmffAngleParams)) {
          if (!fieldOMM) {
            StretchBendContrib *contrib = new StretchBendContrib(
                field, idx[0], idx[1], idx[2], &mmffStbnParams, &mmffAngleParams,
                &mmffBondParams[0], &mmffBondParams[1]);
            field->contribs().push_back(ForceFields::ContribPtr(contrib));
          }
#ifdef RDK_BUILD_WITH_OPENMM
          else {
            fieldOMM->addStretchBendContrib(
                idx[0], idx[1], idx[2], &mmffStbnParams, &mmffAngleParams,
                &mmffBondParams[0], &mmffBondParams[1]);
          }
#endif
          if (mmffMolProperties->getMMFFVerbosity()) {
            unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(idx[0]);
            unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(idx[1]);
            unsigned int kAtomType = mmffMolProperties->getMMFFAtomType(idx[2]);
            const double dist1 = field->distance(idx[0], idx[1]);
            const double dist2 = field->distance(idx[1], idx[2]);
            const RDGeom::Point3D p1((*(points[idx[0]]))[0],
                                     (*(points[idx[0]]))[1],
                                     (*(points[idx[0]]))[2]);
            const RDGeom::Point3D p2((*(points[idx[1]]))[0],
                                     (*(points[idx[1]]))[1],
                                     (*(points[idx[1]]))[2]);
            const RDGeom::Point3D p3((*(points[idx[2]]))[0],
                                     (*(points[idx[2]]))[1],
                                     (*(points[idx[2]]))[2]);
            const double cosTheta =
                MMFF::Utils::calcCosTheta(p1, p2, p3, dist1, dist2);
            const double theta = RAD2DEG * acos(cosTheta);
            const std::pair<double, double> forceConstants =
                MMFF::Utils::calcStbnForceConstants(&mmffStbnParams);
            const std::pair<double, double> stretchBendEnergies =
                MMFF::Utils::calcStretchBendEnergy(
                    dist1 - mmffBondParams[0].r0, dist2 - mmffBondParams[1].r0,
                    theta - mmffAngleParams.theta0, forceConstants);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              if (!isDoubleZero(forceConstants.first)) {
                oStream << std::left << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1 << std::setw(2)
                        << jAtom->getSymbol() << " #" << std::setw(5)
                        << idx[1] + 1 << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1 << std::right
                        << std::setw(5) << iAtomType << std::setw(5)
                        << jAtomType << std::setw(5) << kAtomType
                        << std::setw(6) << stretchBendType << "  " << std::fixed
                        << std::setprecision(3) << std::setw(10) << theta
                        << std::setw(10) << theta - mmffAngleParams.theta0
                        << std::setw(10) << dist1 - mmffBondParams[0].r0
                        << std::setw(10) << dist2 - mmffBondParams[1].r0
                        << std::setw(10) << stretchBendEnergies.first
                        << std::setw(10) << mmffStbnParams.kbaIJK << std::endl;
              }
              if (!isDoubleZero(forceConstants.second)) {
                oStream << std::left << std::setw(2) << kAtom->getSymbol()
                        << " #" << std::setw(5) << idx[2] + 1 << std::setw(2)
                        << jAtom->getSymbol() << " #" << std::setw(5)
                        << idx[1] + 1 << std::setw(2) << iAtom->getSymbol()
                        << " #" << std::setw(5) << idx[0] + 1 << std::right
                        << std::setw(5) << kAtomType << std::setw(5)
                        << jAtomType << std::setw(5) << iAtomType
                        << std::setw(6) << stretchBendType << "  " << std::fixed
                        << std::setprecision(3) << std::setw(10) << theta
                        << std::setw(10) << theta - mmffAngleParams.theta0
                        << std::setw(10) << dist1 - mmffBondParams[0].r0
                        << std::setw(10) << dist2 - mmffBondParams[1].r0
                        << std::setw(10) << stretchBendEnergies.second
                        << std::setw(10) << mmffStbnParams.kbaKJI << std::endl;
              }
            }
            totalStretchBendEnergy +=
                (stretchBendEnergies.first + stretchBendEnergies.second);
          }
        }
        ++j;
      }
      ++i;
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL STRETCH-BEND ENERGY      =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalStretchBendEnergy
            << std::endl;
  }
}

void addOop(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field) {
  addOop(mol, mmffMolProperties, field, 0);
}

void addOop(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field, int ffOpts) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ffOpts);
  unsigned int idx[4];
  unsigned int atomType[4];
  unsigned int n[4];
  const Atom *atom[4];
  ROMol::ADJ_ITER nbrIdx;
  ROMol::ADJ_ITER endNbrs;

  double totalOopBendEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "O U T - O F - P L A N E   B E N D I N G\n\n"
                 "--------------ATOMS---------------         ATOM TYPES        "
                 " OOP                FORCE\n"
                 "  I        J        K        L          I    J    K    L     "
                 "ANGLE    ENERGY   CONSTANT\n"
                 "-------------------------------------------------------------"
                 "-----------------------------" << std::endl;
    }
    points = field->positions();
  }
  for (idx[1] = 0; idx[1] < mol.getNumAtoms(); ++idx[1]) {
    atom[1] = mol.getAtomWithIdx(idx[1]);
    if (atom[1]->getDegree() != 3) {
      continue;
    }
    atomType[1] = mmffMolProperties->getMMFFAtomType(idx[1]);
    boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(atom[1]);
    unsigned int i = 0;
    for (; nbrIdx != endNbrs; ++nbrIdx) {
      atom[i] = mol[*nbrIdx].get();
      idx[i] = atom[i]->getIdx();
      atomType[i] = mmffMolProperties->getMMFFAtomType(idx[i]);
      if (!i) {
        ++i;
      }
      ++i;
    }
    MMFFOop mmffOopParams;
    // if no parameters could be found, we exclude this term (SURDOX02)
    if (!(mmffMolProperties->getMMFFOopBendParams(mol, idx[0], idx[1], idx[2],
                                                  idx[3], mmffOopParams))) {
      continue;
    }
    for (unsigned int i = 0; i < 3; ++i) {
      n[1] = 1;
      switch (i) {
        case 0:
          n[0] = 0;
          n[2] = 2;
          n[3] = 3;
          break;

        case 1:
          n[0] = 0;
          n[2] = 3;
          n[3] = 2;
          break;

        case 2:
          n[0] = 2;
          n[2] = 3;
          n[3] = 0;
          break;
      }
      if (!fieldOMM) {
        OopBendContrib *contrib = new OopBendContrib(
            field, idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]], &mmffOopParams);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
      }
#ifdef RDK_BUILD_WITH_OPENMM
      else {
        fieldOMM->addOopBendContrib(
            idx[n[0]], idx[n[1]], idx[n[2]], idx[n[3]], &mmffOopParams);
      }
#endif
      if (mmffMolProperties->getMMFFVerbosity()) {
        const RDGeom::Point3D p1((*(points[idx[n[0]]]))[0],
                                 (*(points[idx[n[0]]]))[1],
                                 (*(points[idx[n[0]]]))[2]);
        const RDGeom::Point3D p2((*(points[idx[n[1]]]))[0],
                                 (*(points[idx[n[1]]]))[1],
                                 (*(points[idx[n[1]]]))[2]);
        const RDGeom::Point3D p3((*(points[idx[n[2]]]))[0],
                                 (*(points[idx[n[2]]]))[1],
                                 (*(points[idx[n[2]]]))[2]);
        const RDGeom::Point3D p4((*(points[idx[n[3]]]))[0],
                                 (*(points[idx[n[3]]]))[1],
                                 (*(points[idx[n[3]]]))[2]);
        const double chi = MMFF::Utils::calcOopChi(p1, p2, p3, p4);
        const double oopBendEnergy =
            MMFF::Utils::calcOopBendEnergy(chi, mmffOopParams.koop);
        if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
          oStream << std::left << std::setw(2) << atom[0]->getSymbol() << " #"
                  << std::setw(5) << idx[n[0]] + 1 << std::setw(2)
                  << atom[1]->getSymbol() << " #" << std::setw(5)
                  << idx[n[1]] + 1 << std::setw(2) << atom[2]->getSymbol()
                  << " #" << std::setw(5) << idx[n[2]] + 1 << std::setw(2)
                  << atom[3]->getSymbol() << " #" << std::setw(5)
                  << idx[n[3]] + 1 << std::right << std::setw(5)
                  << atomType[n[0]] << std::setw(5) << atomType[n[1]]
                  << std::setw(5) << atomType[n[2]] << std::setw(5)
                  << atomType[n[3]] << std::fixed << std::setprecision(3)
                  << std::setw(10) << chi << std::setw(10) << oopBendEnergy
                  << std::setw(10) << mmffOopParams.koop << std::endl;
        }
        totalOopBendEnergy += oopBendEnergy;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL OUT-OF-PLANE BEND ENERGY =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalOopBendEnergy
            << std::endl;
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
const std::string DefaultTorsionBondSmarts::ds_string = "[!$(*#*)&!D1]~[!$(*#*)&!D1]";
boost::scoped_ptr<const ROMol> DefaultTorsionBondSmarts::ds_instance;
#ifdef BOOST_THREAD_PROVIDES_ONCE_CXX11
boost::once_flag DefaultTorsionBondSmarts::ds_flag;
#else
boost::once_flag DefaultTorsionBondSmarts::ds_flag = BOOST_ONCE_INIT;
#endif

void DefaultTorsionBondSmarts::create() {
  ds_instance.reset(SmartsToMol(ds_string));
}

const ROMol *DefaultTorsionBondSmarts::query() {
#ifdef RDK_THREADSAFE_SSS
  boost::call_once(&create, ds_flag);
#else
  static bool created = false;
  if (!created) {
    created = true;
    create();
  }
#endif
  return ds_instance.get();
}

void addTorsions(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                 ForceFields::ForceField *field,
                 const std::string &torsionBondSmarts) {
  addTorsions(mol, mmffMolProperties, field, 0, torsionBondSmarts);
}

void addTorsions(const ROMol &mol, MMFFMolProperties *mmffMolProperties,
                 ForceFields::ForceField *field, int ffOpts,
                 const std::string &torsionBondSmarts) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);
  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ffOpts);
  ROMol::ADJ_ITER nbr1Idx;
  ROMol::ADJ_ITER end1Nbrs;
  ROMol::ADJ_ITER nbr2Idx;
  ROMol::ADJ_ITER end2Nbrs;
  double totalTorsionEnergy = 0.0;
  RDGeom::PointPtrVect points;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "T O R S I O N A L\n\n"
                 "--------------ATOMS---------------      ---ATOM TYPES---     "
                 "FF     TORSION              -----FORCE Params-----\n"
                 "  I        J        K        L          I    J    K    L   "
                 "CLASS     ANGLE    ENERGY       V1        V2        V3\n"
                 "-------------------------------------------------------------"
                 "-----------------------------------------------------"
              << std::endl;
    }
    points = field->positions();
  }
  std::vector<MatchVectType> matchVect;
  const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
  const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
                       ? defaultQuery : SmartsToMol(torsionBondSmarts);
  TEST_ASSERT(query);
  unsigned int nHits = SubstructMatch(mol, *query, matchVect);
  if (query != defaultQuery) delete query;

  for (unsigned int i = 0; i < nHits; ++i) {
    MatchVectType match = matchVect[i];
    TEST_ASSERT(match.size() == 2);
    int idx2 = match[0].second;
    int idx3 = match[1].second;
    const Bond *bond = mol.getBondBetweenAtoms(idx2, idx3);
    TEST_ASSERT(bond);
    const Atom *jAtom = mol.getAtomWithIdx(idx2);
    const Atom *kAtom = mol.getAtomWithIdx(idx3);
    if (((jAtom->getHybridization() == Atom::SP2) ||
         (jAtom->getHybridization() == Atom::SP3)) &&
        ((kAtom->getHybridization() == Atom::SP2) ||
         (kAtom->getHybridization() == Atom::SP3))) {
      ROMol::OEDGE_ITER beg1, end1;
      boost::tie(beg1, end1) = mol.getAtomBonds(jAtom);
      while (beg1 != end1) {
        const Bond *tBond1 = mol[*beg1].get();
        if (tBond1 != bond) {
          int idx1 = tBond1->getOtherAtomIdx(idx2);
          ROMol::OEDGE_ITER beg2, end2;
          boost::tie(beg2, end2) = mol.getAtomBonds(kAtom);
          while (beg2 != end2) {
            const Bond *tBond2 = mol[*beg2].get();
            if ((tBond2 != bond) && (tBond2 != tBond1)) {
              int idx4 = tBond2->getOtherAtomIdx(idx3);
              // make sure this isn't a three-membered ring:
              if (idx4 != idx1) {
                // we now have a torsion involving atoms (bonds):
                //  bIdx - (tBond1) - idx1 - (bond) - idx2 - (tBond2) - eIdx
                unsigned int torType;
                MMFFTor mmffTorParams;
                if (mmffMolProperties->getMMFFTorsionParams(
                        mol, idx1, idx2, idx3, idx4, torType, mmffTorParams)) {
                  if (!fieldOMM) {
                    TorsionAngleContrib *contrib = new TorsionAngleContrib(
                        field, idx1, idx2, idx3, idx4, &mmffTorParams);
                    field->contribs().push_back(ForceFields::ContribPtr(contrib));
                  }
#ifdef RDK_BUILD_WITH_OPENMM
                  else {
                    fieldOMM->addTorsionAngleContrib(
                        idx1, idx2, idx3, idx4, &mmffTorParams);
                  }
#endif
                  if (mmffMolProperties->getMMFFVerbosity()) {
                    const Atom *iAtom = mol.getAtomWithIdx(idx1);
                    const Atom *lAtom = mol.getAtomWithIdx(idx4);
                    unsigned int iAtomType =
                        mmffMolProperties->getMMFFAtomType(idx1);
                    unsigned int jAtomType =
                        mmffMolProperties->getMMFFAtomType(idx2);
                    unsigned int kAtomType =
                        mmffMolProperties->getMMFFAtomType(idx3);
                    unsigned int lAtomType =
                        mmffMolProperties->getMMFFAtomType(idx4);
                    const RDGeom::Point3D p1((*(points[idx1]))[0],
                                             (*(points[idx1]))[1],
                                             (*(points[idx1]))[2]);
                    const RDGeom::Point3D p2((*(points[idx2]))[0],
                                             (*(points[idx2]))[1],
                                             (*(points[idx2]))[2]);
                    const RDGeom::Point3D p3((*(points[idx3]))[0],
                                             (*(points[idx3]))[1],
                                             (*(points[idx3]))[2]);
                    const RDGeom::Point3D p4((*(points[idx4]))[0],
                                             (*(points[idx4]))[1],
                                             (*(points[idx4]))[2]);
                    const double cosPhi =
                        MMFF::Utils::calcTorsionCosPhi(p1, p2, p3, p4);
                    const double torsionEnergy = MMFF::Utils::calcTorsionEnergy(
                        mmffTorParams.V1, mmffTorParams.V2, mmffTorParams.V3,
                        cosPhi);
                    if (mmffMolProperties->getMMFFVerbosity() ==
                        MMFF_VERBOSITY_HIGH) {
                      oStream
                          << std::left << std::setw(2) << iAtom->getSymbol()
                          << " #" << std::setw(5) << idx1 + 1 << std::setw(2)
                          << jAtom->getSymbol() << " #" << std::setw(5)
                          << idx2 + 1 << std::setw(2) << kAtom->getSymbol()
                          << " #" << std::setw(5) << idx3 + 1 << std::setw(2)
                          << lAtom->getSymbol() << " #" << std::setw(5)
                          << idx4 + 1 << std::right << std::setw(5) << iAtomType
                          << std::setw(5) << jAtomType << std::setw(5)
                          << kAtomType << std::setw(5) << lAtomType
                          << std::setw(6) << torType << "  " << std::fixed
                          << std::setprecision(3) << std::setw(10)
                          << RAD2DEG * acos(cosPhi) << std::setw(10)
                          << torsionEnergy << std::setw(10) << mmffTorParams.V1
                          << std::setw(10) << mmffTorParams.V2 << std::setw(10)
                          << mmffTorParams.V3 << std::endl;
                    }
                    totalTorsionEnergy += torsionEnergy;
                  }
                }
              }
            }
            beg2++;
          }
        }
        beg1++;
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL TORSIONAL ENERGY         =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalTorsionEnergy
            << std::endl;
  }
}

#ifdef RDK_BUILD_WITH_OPENMM
void addNonbonded(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field,
            boost::shared_array<boost::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  checkFFPreconditions(field, mmffMolProperties, ForceFields::USE_OPENMM);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  OpenMMForceField *fieldOMM = getOpenMMForceField(field, ForceFields::USE_OPENMM);
  INT_VECT fragMapping;
  if (ignoreInterfragInteractions) {
    std::vector<ROMOL_SPTR> molFrags =
        MolOps::getMolFrags(mol, true, &fragMapping);
  }
  unsigned int nAtoms = mol.getNumAtoms();
  double totalVdWEnergy = 0.0;
  double totalEleEnergy = 0.0;
  double dielConst = mmffMolProperties->getMMFFDielectricConstant();
  boost::uint8_t dielModel = mmffMolProperties->getMMFFDielectricModel();
  if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
    oStream << "\n"
               "V A N   D E R   W A A L S\n\n"
               "------ATOMS------   ATOM TYPES                               "
               "  WELL\n"
               "  I        J          I    J    DISTANCE   ENERGY     R*     "
               " DEPTH\n"
               "-------------------------------------------------------------"
               "-------" << std::endl;
  }
  const Conformer &conf = mol.getConformer(confId);
  MMFFVdWCollection *mmffVdW = MMFFVdWCollection::getMMFFVdW();
  std::vector<std::pair<int, int> > excl;
  bool haveNonbondedContrib = false;
  for (unsigned int i = 0; i < nAtoms; ++i) {
    const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
    const MMFFVdW *mmffVdWParams = (*mmffVdW)(iAtomType);
    haveNonbondedContrib = true;
    fieldOMM->addNonbondedContrib(i, mmffMolProperties->getMMFFEleTerm()
      ? mmffMolProperties->getMMFFPartialCharge(i) : 0.0,
      mmffMolProperties->getMMFFVdWTerm() ? mmffVdWParams : NULL);
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
        excl.push_back(std::make_pair(i, j));
        continue;
      }
      if (mmffMolProperties->getMMFFVerbosity()
        && getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >= RELATION_1_4) {
        double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        if (dist > nonBondedThresh) {
          continue;
        }
        MMFFVdWRijstarEps mmffVdWConstants;
        if (mmffMolProperties->getMMFFVdWParams(i, j, mmffVdWConstants)) {
          const double vdWEnergy = MMFF::Utils::calcVdWEnergy(
              dist, mmffVdWConstants.R_ij_star, mmffVdWConstants.epsilon);
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
            const Atom *iAtom = mol.getAtomWithIdx(i);
            const Atom *jAtom = mol.getAtomWithIdx(j);
            oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                    << std::setw(5) << i + 1 << std::setw(2)
                    << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType << std::setw(5)
                    << jAtomType << "  " << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist << std::setw(10) << vdWEnergy
                    << std::setw(9) << mmffVdWConstants.R_ij_star << std::setw(9)
                    << mmffVdWConstants.epsilon << std::endl;
          }
          totalVdWEnergy += vdWEnergy;
        }
      }
    }
#if 0
    std::cerr << "particle " << i << ", excl = ";
    for (std::vector<int>::const_iterator it = excl.begin(); it != excl.end(); ++it)
      std::cerr << *it << ((it == excl.end() - 1) ? "\n" : ", ");
#endif
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL VAN DER WAALS ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalVdWEnergy
            << std::endl;
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "E L E C T R O S T A T I C\n\n"
                 "------ATOMS------   ATOM TYPES\n"
                 "  I        J          I    J    DISTANCE   ENERGY\n"
                 "--------------------------------------------------"
              << std::endl;
    }
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (!(mmffMolProperties->getMMFFEleTerm()
        && !isDoubleZero(mmffMolProperties->getMMFFPartialCharge(i))))
        continue;
      for (unsigned int j = i + 1; j < nAtoms; ++j) {
        if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
          continue;
        }
        boost::uint8_t cell = getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j));
        bool is1_4 = (cell == RELATION_1_4);
        if (cell >= RELATION_1_4) {
          double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
          if (dist > nonBondedThresh) {
            continue;
          }
          double chargeTerm = mmffMolProperties->getMMFFPartialCharge(i) *
                              mmffMolProperties->getMMFFPartialCharge(j) /
                              dielConst;
          const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
          const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
          const Atom *iAtom = mol.getAtomWithIdx(i);
          const Atom *jAtom = mol.getAtomWithIdx(j);
          const double eleEnergy = MMFF::Utils::calcEleEnergy(
            i, j, dist, chargeTerm, dielModel, is1_4);
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                    << std::setw(5) << i + 1 << std::setw(2)
                    << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType << std::setw(5)
                    << jAtomType << "  " << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist << std::setw(10) << eleEnergy
                    << std::endl;
          }
          totalEleEnergy += eleEnergy;
        }
      }
    }
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL ELECTROSTATIC ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalEleEnergy
            << std::endl;
  }
  if (haveNonbondedContrib) {
    std::vector<std::pair<int, int> > bonds;
    for (ROMol::ConstBondIterator it = mol.beginBonds(); it != mol.endBonds(); ++it)
      bonds.push_back(std::make_pair((*it)->getBeginAtomIdx(), (*it)->getEndAtomIdx()));
    fieldOMM->addNonbondedExclusionsAndExceptions(excl, bonds);
  }
}

#endif

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addVdW(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field,
            boost::shared_array<boost::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  addVdW(mol, confId, mmffMolProperties, field, 0, neighborMatrix,
         nonBondedThresh, ignoreInterfragInteractions);
}

void addVdW(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field, int ffOpts,
            boost::shared_array<boost::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  INT_VECT fragMapping;
  if (ignoreInterfragInteractions) {
    std::vector<ROMOL_SPTR> molFrags =
        MolOps::getMolFrags(mol, true, &fragMapping);
  }

  unsigned int nAtoms = mol.getNumAtoms();
  double totalVdWEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "V A N   D E R   W A A L S\n\n"
                 "------ATOMS------   ATOM TYPES                               "
                 "  WELL\n"
                 "  I        J          I    J    DISTANCE   ENERGY     R*     "
                 " DEPTH\n"
                 "-------------------------------------------------------------"
                 "-------" << std::endl;
    }
  }
  const Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
        continue;
      }
      if (getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j)) >= RELATION_1_4) {
        double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        if (dist > nonBondedThresh) {
          continue;
        }
        MMFFVdWRijstarEps mmffVdWConstants;
        if (mmffMolProperties->getMMFFVdWParams(i, j, mmffVdWConstants)) {
          VdWContrib *contrib = new VdWContrib(field, i, j, &mmffVdWConstants);
          field->contribs().push_back(ForceFields::ContribPtr(contrib));
          if (mmffMolProperties->getMMFFVerbosity() && (j > i)) {
            const Atom *iAtom = mol.getAtomWithIdx(i);
            const Atom *jAtom = mol.getAtomWithIdx(j);
            const double vdWEnergy = MMFF::Utils::calcVdWEnergy(
                dist, mmffVdWConstants.R_ij_star, mmffVdWConstants.epsilon);
            if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
              unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
              unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
              oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                      << std::setw(5) << i + 1 << std::setw(2)
                      << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                      << std::right << std::setw(5) << iAtomType << std::setw(5)
                      << jAtomType << "  " << std::fixed << std::setprecision(3)
                      << std::setw(9) << dist << std::setw(10) << vdWEnergy
                      << std::setw(9) << mmffVdWConstants.R_ij_star << std::setw(9)
                      << mmffVdWConstants.epsilon << std::endl;
            }
            totalVdWEnergy += vdWEnergy;
          }
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL VAN DER WAALS ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalVdWEnergy
            << std::endl;
  }
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
void addEle(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field,
            boost::shared_array<boost::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  addEle(mol, confId, mmffMolProperties, field, 0, neighborMatrix,
         nonBondedThresh, ignoreInterfragInteractions);
}

void addEle(const ROMol &mol, int confId, MMFFMolProperties *mmffMolProperties,
            ForceFields::ForceField *field, int ffOpts,
            boost::shared_array<boost::uint8_t> neighborMatrix,
            double nonBondedThresh, bool ignoreInterfragInteractions) {
  checkFFPreconditions(field, mmffMolProperties, ffOpts);

  std::ostream &oStream = mmffMolProperties->getMMFFOStream();
  INT_VECT fragMapping;
  if (ignoreInterfragInteractions) {
    std::vector<ROMOL_SPTR> molFrags =
        MolOps::getMolFrags(mol, true, &fragMapping);
  }
  unsigned int nAtoms = mol.getNumAtoms();
  double totalEleEnergy = 0.0;
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << "\n"
                 "E L E C T R O S T A T I C\n\n"
                 "------ATOMS------   ATOM TYPES\n"
                 "  I        J          I    J    DISTANCE   ENERGY\n"
                 "--------------------------------------------------"
              << std::endl;
    }
  }
  const Conformer &conf = mol.getConformer(confId);
  double dielConst = mmffMolProperties->getMMFFDielectricConstant();
  boost::uint8_t dielModel = mmffMolProperties->getMMFFDielectricModel();
  for (unsigned int i = 0; i < nAtoms; ++i) {
    for (unsigned int j = i + 1; j < nAtoms; ++j) {
      if (ignoreInterfragInteractions && (fragMapping[i] != fragMapping[j])) {
        continue;
      }
      boost::uint8_t cell = getTwoBitCell(neighborMatrix, twoBitCellPos(nAtoms, i, j));
      bool is1_4 = (cell == RELATION_1_4);
      if (cell >= RELATION_1_4) {
        if (isDoubleZero(mmffMolProperties->getMMFFPartialCharge(i))
          || isDoubleZero(mmffMolProperties->getMMFFPartialCharge(j)))
          continue;
        double dist = (conf.getAtomPos(i) - conf.getAtomPos(j)).length();
        if (dist > nonBondedThresh) {
          continue;
        }
        double chargeTerm = mmffMolProperties->getMMFFPartialCharge(i) *
                            mmffMolProperties->getMMFFPartialCharge(j) /
                            dielConst;
        EleContrib *contrib = new EleContrib(
          field, i, j, chargeTerm, dielModel, is1_4);
        field->contribs().push_back(ForceFields::ContribPtr(contrib));
        if (mmffMolProperties->getMMFFVerbosity() && (j > i)) {
          const unsigned int iAtomType = mmffMolProperties->getMMFFAtomType(i);
          const unsigned int jAtomType = mmffMolProperties->getMMFFAtomType(j);
          const Atom *iAtom = mol.getAtomWithIdx(i);
          const Atom *jAtom = mol.getAtomWithIdx(j);
          const double eleEnergy = MMFF::Utils::calcEleEnergy(
            i, j, dist, chargeTerm, dielModel, is1_4);
          if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
            oStream << std::left << std::setw(2) << iAtom->getSymbol() << " #"
                    << std::setw(5) << i + 1 << std::setw(2)
                    << jAtom->getSymbol() << " #" << std::setw(5) << j + 1
                    << std::right << std::setw(5) << iAtomType << std::setw(5)
                    << jAtomType << "  " << std::fixed << std::setprecision(3)
                    << std::setw(9) << dist << std::setw(10) << eleEnergy
                    << std::endl;
          }
          totalEleEnergy += eleEnergy;
        }
      }
    }
  }
  if (mmffMolProperties->getMMFFVerbosity()) {
    if (mmffMolProperties->getMMFFVerbosity() == MMFF_VERBOSITY_HIGH) {
      oStream << std::endl;
    }
    oStream << "TOTAL ELECTROSTATIC ENERGY     =" << std::right << std::setw(16)
            << std::fixed << std::setprecision(4) << totalEleEnergy
            << std::endl;
  }
}

}  // end of namespace Tools

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
ForceFields::ForceField *constructForceField(ROMol &mol,
  double nonBondedThresh, int confId, bool ignoreInterfragInteractions) {
  MMFFMolProperties mmffMolProperties(mol);
  
  return constructForceField(mol, &mmffMolProperties, 0,
      nonBondedThresh, confId, ignoreInterfragInteractions);
}

// ------------------------------------------------------------------------
//
//
//
// ------------------------------------------------------------------------
ForceFields::ForceField *constructForceField(ROMol &mol,
  MMFFMolProperties *mmffMolProperties, double nonBondedThresh,
  int confId, bool ignoreInterfragInteractions) {

  return constructForceField(mol, mmffMolProperties, 0,
      nonBondedThresh, confId, ignoreInterfragInteractions);
}

ForceFields::ForceField *constructForceField(ROMol &mol,
  MMFFMolProperties *mmffMolProperties, int ffOpts,
  double nonBondedThresh, int confId, bool ignoreInterfragInteractions) {
  PRECONDITION(mmffMolProperties, "bad MMFFMolProperties");
  PRECONDITION(mmffMolProperties->isValid(),
               "missing atom types - invalid force-field");

  ForceFields::ForceField *ff;
  OpenMMForceField *ffOMM = NULL;
#ifdef RDK_BUILD_WITH_OPENMM
  const PeriodicTable *tbl = PeriodicTable::getTable();
  if (ffOpts & ForceFields::USE_OPENMM) {
    ff = new OpenMMForceField();
    ffOMM = static_cast<OpenMMForceField *>(ff);
  }
  std::vector<OpenMM::Vec3> positionsOMM;
#endif
  if (!ffOMM)
    ff = new ForceFields::ForceField();
  RDKit::MMFF::Tools::checkFFPreconditions(ff, mmffMolProperties, ffOpts);
  // add the atomic positions:
  Conformer &conf = mol.getConformer(confId);
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    RDGeom::Point3D *posPtr = &(conf.getAtomPos(i));
    ff->positions().push_back(posPtr);
#ifdef RDK_BUILD_WITH_OPENMM
    if (ffOMM) {
      positionsOMM.push_back(OpenMM::Vec3(
        posPtr->x * OpenMM::NmPerAngstrom,
        posPtr->y * OpenMM::NmPerAngstrom,
        posPtr->z * OpenMM::NmPerAngstrom));
      ffOMM->getSystem()->addParticle(tbl->getAtomicWeight(
        mol.getAtomWithIdx(i)->getAtomicNum()));
    }
#endif
  }
  ff->initialize();
  if (mmffMolProperties->getMMFFBondTerm()) {
    Tools::addBonds(mol, mmffMolProperties, ff, ffOpts);
  }
  if (mmffMolProperties->getMMFFAngleTerm()) {
    Tools::addAngles(mol, mmffMolProperties, ff, ffOpts);
  }
  if (mmffMolProperties->getMMFFStretchBendTerm()) {
    Tools::addStretchBend(mol, mmffMolProperties, ff, ffOpts);
  }
  if (mmffMolProperties->getMMFFOopTerm()) {
    Tools::addOop(mol, mmffMolProperties, ff, ffOpts);
  }
  if (mmffMolProperties->getMMFFTorsionTerm()) {
    Tools::addTorsions(mol, mmffMolProperties, ff, ffOpts);
  }
  if (mmffMolProperties->getMMFFVdWTerm() ||
      mmffMolProperties->getMMFFEleTerm()) {
#ifdef RDK_BUILD_WITH_OPENMM
    bool haveSeparateVdwEle = !(ffOpts & ForceFields::USE_OPENMM);
    boost::shared_array<boost::uint8_t> neighborMat;
    if (haveSeparateVdwEle || mmffMolProperties->getMMFFVerbosity())
      neighborMat = Tools::buildNeighborMatrix(mol);
    if (!haveSeparateVdwEle)
      Tools::addNonbonded(mol, confId, mmffMolProperties, ff, neighborMat,
                    nonBondedThresh, ignoreInterfragInteractions);
#else
    bool haveSeparateVdwEle = true;
    boost::shared_array<boost::uint8_t> neighborMat =
        Tools::buildNeighborMatrix(mol);
#endif
    if (haveSeparateVdwEle) {
      if (mmffMolProperties->getMMFFVdWTerm()) {
        Tools::addVdW(mol, confId, mmffMolProperties, ff, ffOpts, neighborMat,
                      nonBondedThresh, ignoreInterfragInteractions);
      }
      if (mmffMolProperties->getMMFFEleTerm()) {
        Tools::addEle(mol, confId, mmffMolProperties, ff, ffOpts, neighborMat,
                      nonBondedThresh, ignoreInterfragInteractions);
      }
    }
  }
#ifdef RDK_BUILD_WITH_OPENMM
  if (ffOMM) {
    ffOMM->initializeContext();
    ffOMM->getContext(true)->setPositions(positionsOMM);
  }
#endif

  return ff;
}

#ifdef RDK_BUILD_WITH_OPENMM
OpenMMForceField::OpenMMForceField(bool forceLoadPlugins,
  const std::string &pluginsDir) :
  ForceFields::OpenMMForceField(),
  d_bondStretchForce(NULL),
  d_angleBendForce(NULL),
  d_stretchBendForce(NULL),
  d_torsionAngleForce(NULL),
  d_oopBendForce(NULL),
  d_nonbondedForce(NULL) {
  Tools::OpenMMPlugins::instance()->load(forceLoadPlugins, pluginsDir);
}

void OpenMMForceField::cloneSystemTo(OpenMM::System& other) const {
  ForceFields::OpenMMForceField::cloneSystemTo(other);
  for (int i = 0; i < d_openmmSystem->getNumParticles(); ++i)
    other.addParticle(d_openmmSystem->getParticleMass(i));
  for (int i = 0; i < d_openmmSystem->getNumForces(); ++i) {
    const OpenMM::Force &f = d_openmmSystem->getForce(i);
    const OpenMM::MMFFBondForce *bondForce = dynamic_cast<const OpenMM::MMFFBondForce *>(&f);
    const OpenMM::MMFFAngleForce *angleForce = dynamic_cast<const OpenMM::MMFFAngleForce *>(&f);
    const OpenMM::MMFFStretchBendForce *stbnForce = dynamic_cast<const OpenMM::MMFFStretchBendForce *>(&f);
    const OpenMM::MMFFTorsionForce *torsionForce = dynamic_cast<const OpenMM::MMFFTorsionForce *>(&f);
    const OpenMM::MMFFOutOfPlaneBendForce *oopForce = dynamic_cast<const OpenMM::MMFFOutOfPlaneBendForce *>(&f);
    const OpenMM::MMFFNonbondedForce *nonbondedForce = dynamic_cast<const OpenMM::MMFFNonbondedForce *>(&f);
    if (bondForce) {
      OpenMM::MMFFBondForce *bondForceOther = new OpenMM::MMFFBondForce();
      bondForceOther->setMMFFGlobalBondCubic(bondForce->getMMFFGlobalBondCubic());
      bondForceOther->setMMFFGlobalBondQuartic(bondForce->getMMFFGlobalBondQuartic());
      bondForceOther->setUsesPeriodicBoundaryConditions(bondForce->usesPeriodicBoundaryConditions());
      for (int j = 0; j < bondForce->getNumBonds(); ++j) {
        int particle1;
        int particle2;
        double length;
        double quadraticK;
        bondForce->getBondParameters(j, particle1, particle2, length, quadraticK);
        bondForceOther->addBond(particle1, particle2, length, quadraticK);
      }
      other.addForce(bondForceOther);
    }
    else if (angleForce) {
      OpenMM::MMFFAngleForce *angleForceOther = new OpenMM::MMFFAngleForce();
      angleForceOther->setMMFFGlobalAngleCubic(angleForce->getMMFFGlobalAngleCubic());
      angleForceOther->setUsesPeriodicBoundaryConditions(angleForce->usesPeriodicBoundaryConditions());
      for (int j = 0; j < angleForce->getNumAngles(); ++j) {
        int particle1;
        int particle2;
        int particle3;
        double length;
        double quadraticK;
        bool isLinear;
        angleForce->getAngleParameters(j, particle1, particle2, particle3, length, quadraticK, isLinear);
        angleForceOther->addAngle(particle1, particle2, particle3, length, quadraticK, isLinear);
      }
      other.addForce(angleForceOther);
    }
    else if (stbnForce) {
      OpenMM::MMFFStretchBendForce *stbnForceOther = new OpenMM::MMFFStretchBendForce();
      stbnForceOther->setUsesPeriodicBoundaryConditions(stbnForce->usesPeriodicBoundaryConditions());
      for (int j = 0; j < stbnForce->getNumStretchBends(); ++j) {
        int particle1;
        int particle2;
        int particle3;
        double lengthAB;
        double lengthCB;
        double angle;
        double k1;
        double k2;
        stbnForce->getStretchBendParameters(j, particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
        stbnForceOther->addStretchBend(particle1, particle2, particle3, lengthAB, lengthCB, angle, k1, k2);
      }
      other.addForce(stbnForceOther);
    }
    else if (torsionForce) {
      OpenMM::MMFFTorsionForce *torsionForceOther = new OpenMM::MMFFTorsionForce();
      torsionForceOther->setUsesPeriodicBoundaryConditions(torsionForce->usesPeriodicBoundaryConditions());
      for (int j = 0; j < torsionForce->getNumTorsions(); ++j) {
        int particle1;
        int particle2;
        int particle3;
        int particle4;
        double c1;
        double c2;
        double c3;
        torsionForce->getTorsionParameters(j, particle1, particle2, particle3, particle4, c1, c2, c3);
        torsionForceOther->addTorsion(particle1, particle2, particle3, particle4, c1, c2, c3);
      }
      other.addForce(torsionForceOther);
    }
    else if (oopForce) {
      OpenMM::MMFFOutOfPlaneBendForce *oopForceOther = new OpenMM::MMFFOutOfPlaneBendForce();
      oopForceOther->setUsesPeriodicBoundaryConditions(oopForce->usesPeriodicBoundaryConditions());
      for (int j = 0; j < oopForce->getNumOutOfPlaneBends(); ++j) {
        int particle1;
        int particle2;
        int particle3;
        int particle4;
        double k;
        oopForce->getOutOfPlaneBendParameters(j, particle1, particle2, particle3, particle4, k);
        oopForceOther->addOutOfPlaneBend(particle1, particle2, particle3, particle4, k);
      }
      other.addForce(oopForceOther);
    }
    else if (nonbondedForce) {
      OpenMM::MMFFNonbondedForce *nonbondedForceOther = new OpenMM::MMFFNonbondedForce();
      nonbondedForceOther->setNonbondedMethod(nonbondedForce->getNonbondedMethod());
      nonbondedForceOther->setCutoffDistance(nonbondedForce->getCutoffDistance());
      nonbondedForceOther->setUseSwitchingFunction(nonbondedForce->getUseSwitchingFunction());
      nonbondedForceOther->setSwitchingDistance(nonbondedForce->getSwitchingDistance());
      nonbondedForceOther->setReactionFieldDielectric(nonbondedForce->getReactionFieldDielectric());
      nonbondedForceOther->setEwaldErrorTolerance(nonbondedForce->getEwaldErrorTolerance());
      nonbondedForceOther->setUseDispersionCorrection(nonbondedForce->getUseDispersionCorrection());
      nonbondedForceOther->setReciprocalSpaceForceGroup(nonbondedForce->getReciprocalSpaceForceGroup());
      double alpha;
      int nx;
      int ny;
      int nz;
      nonbondedForce->getPMEParameters(alpha, nx, ny, nz);
      nonbondedForceOther->setPMEParameters(alpha, nx, ny, nz);
      for (int j = 0; j < nonbondedForce->getNumParticles(); ++j) {
        double charge;
        double sigma;
        double G_t_alpha;
        double alpha_d_N;
        char vdwDA;
        nonbondedForce->getParticleParameters(j, charge, sigma, G_t_alpha, alpha_d_N, vdwDA);
        nonbondedForceOther->addParticle(charge, sigma, G_t_alpha, alpha_d_N, vdwDA);
      }
      for (int j = 0; j < nonbondedForce->getNumExceptions(); ++j) {
        int particle1;
        int particle2;
        double chargeProd;
        double sigma;
        double epsilon;
        nonbondedForce->getExceptionParameters(j, particle1, particle2, chargeProd, sigma, epsilon);
        nonbondedForceOther->addException(particle1, particle2, chargeProd, sigma, epsilon);
      }
      other.addForce(nonbondedForceOther);
    }
  }
}

void OpenMMForceField::addBondStretchContrib(unsigned int idx1,
  unsigned int idx2, const MMFFBond *mmffBondParams) {
  static const double c1B = MDYNE_A_TO_KCAL_MOL * OpenMM::KJPerKcal
    * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
  if (!d_bondStretchForce) {
    d_bondStretchForce = MMFF::Utils::getOpenMMBondStretchForce();
    d_openmmSystem->addForce(d_bondStretchForce);
  }
  d_bondStretchForce->addBond(idx1, idx2, mmffBondParams->r0 * OpenMM::NmPerAngstrom, c1B * mmffBondParams->kb);
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_bondStretchForce->addBond(" << idx1 << ", " << idx2 << ", "
    << mmffBondParams->r0 << ", "
    << MDYNE_A_TO_KCAL_MOL * mmffBondParams->kb << ");" << std::endl;
#endif
}

void OpenMMForceField::addAngleBendContrib(unsigned int idx1,
  unsigned int idx2, unsigned int idx3, const MMFFAngle *mmffAngleParams,
  const MMFFProp *mmffPropParamsCentralAtom) {
  static const double c1AB = MDYNE_A_TO_KCAL_MOL * OpenMM::KJPerKcal;
  if (!d_angleBendForce) {
    d_angleBendForce = MMFF::Utils::getOpenMMAngleBendForce();
    d_openmmSystem->addForce(d_angleBendForce);
  }
  d_angleBendForce->addAngle(idx1, idx2, idx3, mmffAngleParams->theta0 * OpenMM::RadiansPerDegree,
    c1AB * mmffAngleParams->ka, MMFF::Utils::isLinear(mmffPropParamsCentralAtom));
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_angleBendForce->addAngle(" << idx1 << ", " << idx2 << ", " << idx3 << ", "
    << mmffAngleParams->theta0 << ", "
    << MDYNE_A_TO_KCAL_MOL * mmffAngleParams->ka << ", "
    << (MMFF::Utils::isLinear(mmffPropParamsCentralAtom) ? "true" : "false") << ");" << std::endl;
#endif
}

void OpenMMForceField::addStretchBendContrib(unsigned int idx1,
  unsigned int idx2, unsigned int idx3, const MMFFStbn *mmffStbnParams,
  const MMFFAngle *mmffAngleParams, const MMFFBond *mmffBondParams1,
  const MMFFBond *mmffBondParams2) {
  static const double c1SB = MDYNE_A_TO_KCAL_MOL
    * OpenMM::AngstromsPerNm * OpenMM::KJPerKcal; 
  if (!d_stretchBendForce) {
    d_stretchBendForce = MMFF::Utils::getOpenMMStretchBendForce();
    d_openmmSystem->addForce(d_stretchBendForce);
  }
  std::pair<double, double> forceConstants =
    MMFF::Utils::calcStbnForceConstants(mmffStbnParams);
  d_stretchBendForce->addStretchBend(idx1, idx2, idx3,
    mmffBondParams1->r0 * OpenMM::NmPerAngstrom,
    mmffBondParams2->r0 * OpenMM::NmPerAngstrom,
    mmffAngleParams->theta0 * OpenMM::RadiansPerDegree,
    c1SB * forceConstants.first, c1SB * forceConstants.second);
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_stretchBendForce->addStretchBend("
    << idx1 << ", " << idx2 << ", " << idx3 << ", "
    << mmffBondParams1->r0 << ", "
    << mmffBondParams2->r0 << ", "
    << mmffAngleParams->theta0 << ", "
    << MDYNE_A_TO_KCAL_MOL * forceConstants.first
    << ", " << MDYNE_A_TO_KCAL_MOL * forceConstants.second << ");" << std::endl;
#endif
}

void OpenMMForceField::addTorsionAngleContrib(unsigned int idx1,
  unsigned int idx2, unsigned int idx3, unsigned int idx4,
  const MMFFTor *mmffTorParams) {
  if (!d_torsionAngleForce) {
    d_torsionAngleForce = MMFF::Utils::getOpenMMTorsionAngleForce();
    d_openmmSystem->addForce(d_torsionAngleForce);
  }
  d_torsionAngleForce->addTorsion(idx1, idx2, idx3, idx4, mmffTorParams->V1 * OpenMM::KJPerKcal,
    mmffTorParams->V2 * OpenMM::KJPerKcal, mmffTorParams->V3 * OpenMM::KJPerKcal);
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_torsionAngleForce->addTorsion(" << idx1 << ", " << idx2 << ", "
    << idx3 << ", " << idx4 << ", " << mmffTorParams->V1 << ", "
    << mmffTorParams->V2 << ", "
    << mmffTorParams->V3 << ");" << std::endl;
#endif
}

void OpenMMForceField::addOopBendContrib(unsigned int idx1,
  unsigned int idx2, unsigned int idx3,
  unsigned int idx4, const MMFFOop *mmffOopParams) {
  static const double c2OOP = MDYNE_A_TO_KCAL_MOL * OpenMM::KJPerKcal;
  if (!d_oopBendForce) {
    d_oopBendForce = MMFF::Utils::getOpenMMOopBendForce();
    d_openmmSystem->addForce(d_oopBendForce);
  }
  d_oopBendForce->addOutOfPlaneBend(idx1, idx4, idx3, idx2,
    c2OOP * mmffOopParams->koop);
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_oopBendForce->addOutOfPlaneBend(" << idx1 << ", "
    << idx4 << ", " << idx3 << ", " << idx2 << ", "
    << MDYNE_A_TO_KCAL_MOL * mmffOopParams->koop
    << ");" << std::endl;
#endif
}

void OpenMMForceField::addNonbondedContrib(unsigned int idx,
  double charge, const MMFFVdW *mmffVdWParams) {
  if (!d_nonbondedForce) {
    d_nonbondedForce = MMFF::Utils::getOpenMMNonbondedForce();
    d_nonbondedForce->setUseDispersionCorrection(false);
    d_openmmSystem->addForce(d_nonbondedForce);
  }
  double sigma;
  double G_t_alpha;
  double alpha_d_N;
  char vdwDA;
  if (mmffVdWParams) {
    sigma = mmffVdWParams->R_star * OpenMM::NmPerAngstrom;
    G_t_alpha = mmffVdWParams->G_i * mmffVdWParams->alpha_i;
    alpha_d_N = mmffVdWParams->alpha_i / mmffVdWParams->N_i;
    vdwDA = static_cast<char>(mmffVdWParams->DA);
  }
  else {
    sigma = 1.0;
    G_t_alpha = 0.0;
    alpha_d_N = 1.0;
    vdwDA = '-';
  }
  int openmmIdx = d_nonbondedForce->addParticle(charge, sigma, G_t_alpha, alpha_d_N, vdwDA);
  if (openmmIdx != static_cast<int>(idx)) {
    std::stringstream ss;
    ss << "RDKit idx (" << idx << ") and OpenMM idx (" << openmmIdx
      <<") for this particle differ";
    throw std::runtime_error(ss.str());
  }
#ifdef PRINT_MMFF_FORCES
  std::cerr << "d_nonbondedForce->addParticle(" << charge << ", "
    << mmffVdWParams->R_star << ", " << G_t_alpha << ", " << alpha_d_N << ", '" << vdwDA
    << "');" << std::endl;
#endif
}

void OpenMMForceField::addNonbondedExclusionsAndExceptions(
  const std::vector<std::pair<int, int> > &exclusions,
  const std::vector<std::pair<int, int> > &bonds) {
  PRECONDITION(d_nonbondedForce, "cannot add exclusions/exceptions if nonbondedForce is NULL");
  d_nonbondedForce->createExceptionsFromBonds(bonds);
  for (std::vector<std::pair<int, int> >::const_iterator it = exclusions.begin();
    it != exclusions.end(); ++it)
    d_nonbondedForce->addException(it->first, it->second, 0.0, 1.0, 0.0);
}

void OpenMMForceField::setCutoffDistance(double distance) {
  const double distanceNm = distance * OpenMM::NmPerAngstrom;
  d_nonbondedForce->setCutoffDistance(distanceNm);
}

void OpenMMForceField::setNonbondedPeriodic(bool periodic) {
  OpenMM::MMFFNonbondedForce::NonbondedMethod nonbondedMethod = (periodic
    ? OpenMM::MMFFNonbondedForce::CutoffPeriodic : OpenMM::MMFFNonbondedForce::NoCutoff);
  d_nonbondedForce->setNonbondedMethod(nonbondedMethod);
}

OpenMMForceField *constructOpenMMForceField(ROMol &mol,
  double nonBondedThresh, int confId, bool ignoreInterfragInteractions) {
  MMFFMolProperties mmffMolProperties(mol);

  return constructOpenMMForceField(mol, &mmffMolProperties,
    nonBondedThresh, confId, ignoreInterfragInteractions);
}

OpenMMForceField *constructOpenMMForceField(ROMol &mol,
  MMFFMolProperties *mmffMolProperties, double nonBondedThresh, int confId,
  bool ignoreInterfragInteractions) {

  return static_cast<OpenMMForceField *>(
    constructForceField(mol, mmffMolProperties, ForceFields::USE_OPENMM,
    nonBondedThresh, confId, ignoreInterfragInteractions));
}
#endif

}
}
