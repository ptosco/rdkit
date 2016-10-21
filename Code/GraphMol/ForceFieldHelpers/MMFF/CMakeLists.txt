rdkit_test(testMMFFForceFieldHelpers testMMFFHelpers.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse
           SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment
           GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS})
if(RDK_BUILD_WITH_OPENMM)
rdkit_test(testMMFFForceFieldHelpers_OpenMM testMMFFHelpers_OpenMM.cpp
           LINK_LIBRARIES 
           ForceFieldHelpers
           DistGeomHelpers DistGeometry FileParsers MolTransforms SmilesParse
           SubstructMatch ForceField MolAlign Optimizer EigenSolvers Alignment
           GraphMol RDGeometryLib RDGeneral ${RDKit_THREAD_LIBS} ${OPENMM_LIBRARIES})
endif(RDK_BUILD_WITH_OPENMM)
