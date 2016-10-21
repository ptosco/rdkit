// $Id$
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
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/utils.h>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <ForceField/ForceField.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <boost/math/special_functions/round.hpp>

using namespace RDKit;

void testMMFFBuilder1() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder tools." << std::endl;

  ROMol *mol;

  boost::shared_array<boost::uint8_t> nbrMat;

  mol = SmilesToMol("Cc1c(N)c(O)c(F)c(Cl)c1");
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol) >= 0);
  MMFF::MMFFMolProperties *mmffMolProperties =
      new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  {
    ForceFields::ForceField *field = MMFF::constructForceField(*mol, mmffMolProperties, ForceFields::USE_RDK);
    TEST_ASSERT(field);
    std::cerr << "RDKit MMFF energy = " << field->calcEnergy() << std::endl;
    delete field;
  }
  {
    MMFF::OpenMMForceField *field = static_cast<MMFF::OpenMMForceField *>(
      MMFF::constructForceField(*mol, mmffMolProperties, ForceFields::USE_OPENMM));
    TEST_ASSERT(field);
    std::cerr << "OpenMM MMFF energy = " << field->calcEnergy() << std::endl;
    delete field;
  }
  delete mol;
  delete mmffMolProperties;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testMMFFBuilder2() {
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Testing MMFF builder tools." << std::endl;

  ROMol *mol;

  boost::shared_array<boost::uint8_t> nbrMat;

  mol = SmilesToMol("CCCC");
  TEST_ASSERT(DGeomHelpers::EmbedMolecule(*mol) >= 0);
  MMFF::MMFFMolProperties *mmffMolProperties =
      new MMFF::MMFFMolProperties(*mol);
  TEST_ASSERT(mmffMolProperties);
  TEST_ASSERT(mmffMolProperties->isValid());
  mmffMolProperties->setMMFFVerbosity(MMFF::MMFF_VERBOSITY_HIGH);
  {
    ForceFields::ForceField *field = MMFF::constructForceField(*mol, mmffMolProperties, ForceFields::USE_RDK);
    TEST_ASSERT(field);
    std::cerr << "RDKit MMFF energy = " << field->calcEnergy() << std::endl;
    delete field;
  }
  {
    MMFF::OpenMMForceField *field = static_cast<MMFF::OpenMMForceField *>(
      MMFF::constructForceField(*mol, mmffMolProperties, ForceFields::USE_OPENMM));
    TEST_ASSERT(field);
    std::cerr << "OpenMM MMFF energy = " << field->calcEnergy() << std::endl;
    delete field;
  }
  delete mol;
  delete mmffMolProperties;
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//
//-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
int main() {
  RDLog::InitLogs();
  testMMFFBuilder2();
}
