// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <ForceField/ForceField.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>
#include <GraphMol/ForceFieldHelpers/MMFF/Builder.h>
#include <ForceField/MMFF/Params.h>
#include <GraphMol/Trajectory/Snapshot.h>
#include <boost/python/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <algorithm>
#include <Geometry/point.h>

namespace python = boost::python;

namespace ForceFields {
class PyForceField {
 public:
  PyForceField(ForceField *f) : field(f){};

  ~PyForceField() {
    // std::cerr << " *** destroy PyForce field " << std::endl;
    field.reset();
    // std::cerr << " ***       reset DONE" << std::endl;
    extraPoints.clear();
    // std::cerr << " *** destroy PyForce field DONE" << std::endl;
  }

  int addExtraPoint(double x, double y, double z, bool fixed = true) {
    RDGeom::Point3D *pt = new RDGeom::Point3D(x, y, z);
    PRECONDITION(this->field, "no force field");
    this->extraPoints.push_back(boost::shared_ptr<RDGeom::Point3D>(pt));
    unsigned int ptIdx = this->extraPoints.size() - 1;
    RDGeom::Point3D *ptr = this->extraPoints[ptIdx].get();
    this->field->positions().push_back(ptr);
    int idx = this->field->positions().size();
    if (fixed) {
      this->field->fixedPoints().push_back(idx - 1);
    }
    return idx;
  }

  double calcEnergy() {
    PRECONDITION(this->field, "no force field");
    return this->field->calcEnergy();
  }

  int minimize(int maxIts, double forceTol, double energyTol) {
    PRECONDITION(this->field, "no force field");
    return this->field->minimize(maxIts, forceTol, energyTol);
  }

  python::tuple minimizeTrajectory(unsigned int snapshotFreq, int maxIts, double forceTol, double energyTol);

  void initialize() {
    PRECONDITION(this->field, "no force field");
    this->field->initialize();
  }

  // private:
  std::vector<boost::shared_ptr<RDGeom::Point3D> > extraPoints;
  boost::shared_ptr<ForceField> field;
};

#ifdef RDK_BUILD_WITH_OPENMM
class PyOpenMMForceField {
  public:
    PyOpenMMForceField(OpenMMForceField *f) :
      fieldOMM(f) {};

    ~PyOpenMMForceField() {
      // std::cerr << " *** destroy PyOpenMMForce field " << std::endl;
      fieldOMM.reset();
    }

    OpenMM::System *getSystem() const {
      checkFieldOMM();
      return fieldOMM->getSystem();
    }

    OpenMM::Context *getContext(bool throwIfNull = false) const {
      checkFieldOMM();
      return fieldOMM->getContext(throwIfNull);
    }

    OpenMM::Integrator *getIntegrator() const {
      checkFieldOMM();
      return fieldOMM->getIntegrator();
    }
    
    void setIntegrator(OpenMM::Integrator *integrator) const {
      checkFieldOMM();
      fieldOMM->setIntegrator(integrator);
    }
    
    void initialize() {
      checkFieldOMM();
      fieldOMM->initialize();
    }

    void initializeContext(std::string platformName = std::string(),
      python::dict properties = python::dict()) {
      checkFieldOMM();
      if (!platformName.size())
        fieldOMM->initializeContext();
      else {
        python::list keys = properties.keys();
        std::map<std::string, std::string> propertiesMap;
        for (unsigned int i = 0; i < len(keys); ++i) {
          python::object value = properties[keys[i]];
          if (value)
            propertiesMap[python::extract<std::string>(keys[i])] =
              python::extract<std::string>(value);
        }
        fieldOMM->initializeContext(platformName, propertiesMap);
      }
    }

    double calcEnergy() const {
      checkFieldOMM();
      return fieldOMM->calcEnergy();
    }

    int minimize(int maxIts, double forceTol, double energyTol) {
      checkFieldOMM();
      return fieldOMM->minimize(maxIts, forceTol, energyTol);
    }

    // private:
    boost::shared_ptr<OpenMMForceField> fieldOMM;
  private:
    void checkFieldOMM() const {
      PRECONDITION(fieldOMM.get(), "no force field");
    }
};

class PyMMFFOpenMMForceField : public PyOpenMMForceField {
  public:
    PyMMFFOpenMMForceField(RDKit::MMFF::OpenMMForceField *f) :
      PyOpenMMForceField(f) {};
    python::list loadedPlugins() {
      python::list lp;
      RDKit::MMFF::OpenMMForceField *f = upcast();
      for (std::vector<std::string>::const_iterator it = f->loadedPlugins().begin();
        it != f->loadedPlugins().end() ; ++it)
        lp.append(*it);
      return lp;
    }
    python::list failedPlugins() {
      python::list fp;
      RDKit::MMFF::OpenMMForceField *f = upcast();
      for (std::vector<std::string>::const_iterator it = f->failedPlugins().begin();
        it != f->failedPlugins().end() ; ++it)
        fp.append(*it);
      return fp;
    }
  private:
    RDKit::MMFF::OpenMMForceField *upcast() const {
      PRECONDITION(fieldOMM.get(), "no force field");
      RDKit::MMFF::OpenMMForceField *super = dynamic_cast<RDKit::MMFF::OpenMMForceField *>(fieldOMM.get());
      PRECONDITION(super, "not a MMFF OpenMM force field");
      return super;
    }
};
#endif

class PyMMFFMolProperties {
 public:
  PyMMFFMolProperties(RDKit::MMFF::MMFFMolProperties *mp)
      : mmffMolProperties(mp){};
  ~PyMMFFMolProperties(){};

  unsigned int getMMFFAtomType(unsigned int idx) {
    return (unsigned int)(mmffMolProperties->getMMFFAtomType(idx));
  };
  double getMMFFFormalCharge(unsigned int idx) {
    return mmffMolProperties->getMMFFFormalCharge(idx);
  };
  double getMMFFPartialCharge(unsigned int idx) {
    return mmffMolProperties->getMMFFPartialCharge(idx);
  };
  PyObject *getMMFFBondStretchParams(const RDKit::ROMol &mol,
                                     const unsigned int idx1,
                                     const unsigned int idx2);
  PyObject *getMMFFAngleBendParams(const RDKit::ROMol &mol,
                                   const unsigned int idx1,
                                   const unsigned int idx2,
                                   const unsigned int idx3);
  PyObject *getMMFFStretchBendParams(const RDKit::ROMol &mol,
                                     const unsigned int idx1,
                                     const unsigned int idx2,
                                     const unsigned int idx3);
  PyObject *getMMFFTorsionParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4);
  PyObject *getMMFFOopBendParams(const RDKit::ROMol &mol,
                                 const unsigned int idx1,
                                 const unsigned int idx2,
                                 const unsigned int idx3,
                                 const unsigned int idx4);
  PyObject *getMMFFVdWParams(const unsigned int idx1, const unsigned int idx2);
  void setMMFFDielectricModel(boost::uint8_t dielModel) {
    mmffMolProperties->setMMFFDielectricModel(dielModel);
  };
  void setMMFFDielectricConstant(double dielConst) {
    mmffMolProperties->setMMFFDielectricConstant(dielConst);
  };
  void setMMFFBondTerm(bool state) {
    mmffMolProperties->setMMFFBondTerm(state);
  };
  void setMMFFAngleTerm(const bool state) {
    mmffMolProperties->setMMFFAngleTerm(state);
  };
  void setMMFFStretchBendTerm(const bool state) {
    mmffMolProperties->setMMFFStretchBendTerm(state);
  };
  void setMMFFOopTerm(const bool state) {
    mmffMolProperties->setMMFFOopTerm(state);
  };
  void setMMFFTorsionTerm(const bool state) {
    mmffMolProperties->setMMFFTorsionTerm(state);
  };
  void setMMFFVdWTerm(const bool state) {
    mmffMolProperties->setMMFFVdWTerm(state);
  };
  void setMMFFEleTerm(const bool state) {
    mmffMolProperties->setMMFFEleTerm(state);
  };
  void setMMFFVariant(const std::string &mmffVariant) {
    mmffMolProperties->setMMFFVariant(mmffVariant);
  };
  void setMMFFVerbosity(unsigned int verbosity) {
    mmffMolProperties->setMMFFVerbosity(verbosity);
  };
  boost::shared_ptr<RDKit::MMFF::MMFFMolProperties> mmffMolProperties;
};
PyObject *getUFFBondStretchParams(const RDKit::ROMol &mol,
                                  const unsigned int idx1,
                                  const unsigned int idx2);
PyObject *getUFFAngleBendParams(const RDKit::ROMol &mol,
                                const unsigned int idx1,
                                const unsigned int idx2,
                                const unsigned int idx3);
PyObject *getUFFTorsionParams(const RDKit::ROMol &mol, const unsigned int idx1,
                              const unsigned int idx2, const unsigned int idx3,
                              const unsigned int idx4);
PyObject *getUFFInversionParams(const RDKit::ROMol &mol,
                                const unsigned int idx1,
                                const unsigned int idx2,
                                const unsigned int idx3,
                                const unsigned int idx4);
PyObject *getUFFVdWParams(const RDKit::ROMol &mol, const unsigned int idx1,
                          const unsigned int idx2);
}  // end namespace ForceFields
