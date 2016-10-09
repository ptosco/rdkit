// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "BondStretch.h"
#include "Params.h"
#include <sstream>
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#ifdef RDK_BUILD_WITH_OPENMM
#include <OpenMM.h>
#include <openmm/Units.h>
#endif

namespace ForceFields {
namespace MMFF {
namespace Utils {

static const double c1 = 0.5 * MDYNE_A_TO_KCAL_MOL;
static const double cs = -2.0;
static const double c3 = 7.0 / 12.0 * cs * cs;

double calcBondRestLength(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->r0;
}

double calcBondForceConstant(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->kb;
}

double calcBondStretchEnergy(const double r0, const double kb,
                             const double distance) {
  double distTerm = distance - r0;
  double distTerm2 = distTerm * distTerm;

  return (c1 * kb * distTerm2 *
         (1.0 + cs * distTerm + c3 * distTerm2));
}

#ifdef RDK_BUILD_WITH_OPENMM
static const double c1OMM = c1 * OpenMM::KJPerKcal
  * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
static const double csOMM = cs * OpenMM::AngstromsPerNm;
static const double c3OMM = c3 * OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;

OpenMM::CustomBondForce *getOpenMMBondStretchForce() {
  std::stringstream bf;
  bf << c1OMM << "*kb*(r-r0)^2*(1.0+"
     << csOMM << "*(r-r0)+" << c3OMM << "*(r-r0)^2)";
  OpenMM::CustomBondForce *res = new OpenMM::CustomBondForce(bf.str());
  res->addPerBondParameter("kb");
  res->addPerBondParameter("r0");
  
  return res;
}
#endif
}  // end of namespace Utils

BondStretchContrib::BondStretchContrib(ForceField *owner,
                                       const unsigned int idx1,
                                       const unsigned int idx2,
                                       const MMFFBond *mmffBondParams) {
  PRECONDITION(owner, "bad owner");
  URANGE_CHECK(idx1, owner->positions().size() - 1);
  URANGE_CHECK(idx2, owner->positions().size() - 1);

  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_r0 = mmffBondParams->r0;
  d_kb = mmffBondParams->kb;
}

double BondStretchContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  return Utils::calcBondStretchEnergy(
      d_r0, d_kb, dp_forceField->distance(d_at1Idx, d_at2Idx, pos));
}

void BondStretchContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);

  double *at1Coords = &(pos[3 * d_at1Idx]);
  double *at2Coords = &(pos[3 * d_at2Idx]);
  double *g1 = &(grad[3 * d_at1Idx]);
  double *g2 = &(grad[3 * d_at2Idx]);
  double const cs = -2.0;
  double const c1 = MDYNE_A_TO_KCAL_MOL;
  double const c3 = 7.0 / 12.0;
  double distTerm = dist - d_r0;
  double dE_dr =
      c1 * d_kb * distTerm *
      (1.0 + 1.5 * cs * distTerm + 2.0 * c3 * cs * cs * distTerm * distTerm);
  double dGrad;
  for (unsigned int i = 0; i < 3; ++i) {
    dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
                          : d_kb * 0.01);
    g1[i] += dGrad;
    g2[i] -= dGrad;
  }
}
}
}
