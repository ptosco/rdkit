if(RDK_BUILD_PYTHON_STUBS)
  string(CONCAT RUN_GEN_RDKIT_STUBS_PY "execute_process("
      "COMMAND "
      "${PYTHON_EXECUTABLE} \"${CMAKE_SOURCE_DIR}/Scripts/gen_rdkit_stubs.py\" "
      "\"${CMAKE_SOURCE_DIR}\" \"${PYTHON_INSTDIR}\" "
      "WORKING_DIRECTORY \"${CMAKE_SOURCE_DIR}\" "
      "OUTPUT_FILE \"${CMAKE_CURRENT_BINARY_DIR}/gen_rdkit_stubs.out\" "
      "ERROR_FILE \"${CMAKE_CURRENT_BINARY_DIR}/gen_rdkit_stubs.err\""
    ")"
  )
  install(CODE "${RUN_GEN_RDKIT_STUBS_PY}" COMPONENT python)
else(RDK_BUILD_PYTHON_STUBS)
  install(DIRECTORY "${CMAKE_SOURCE_DIR}/rdkit-stubs"
    DESTINATION ${PYTHON_INSTDIR}
    COMPONENT python
    PATTERN "CMakeLists.txt" EXCLUDE
  )
endif(RDK_BUILD_PYTHON_STUBS)
