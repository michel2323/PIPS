#add_definitions(-DTIMING -DSTOCH_TESTING) # timing output
include_directories(Core/Abstract Core/Vector Core/Utilities Core/QpSolvers Core/QpGen
  Core/SparseLinearAlgebra Core/DenseLinearAlgebra Core/Readers
  Core/LinearSolvers/Ma27Solver Core/LinearSolvers/Ma57Solver
  Core/LinearSolvers/Ma86Solver Core/LinearSolvers/PardisoSolver Core/LinearSolvers/BiCGStabSolver)
include_directories(Core/StochLinearAlgebra Core/QpStoch)
add_subdirectory(Core)

set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")

if (HAVE_MA57 AND HAVE_METIS)
    add_executable(qpgen-sparse-ma57-mehrotra.exe Core/QpGen/QpGenSparseMa57MehrotraDriver.C)
    target_link_libraries(qpgen-sparse-ma57-mehrotra.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

    # add_executable(pipsipm-unitcommit.exe Drivers/unitCommitment.cpp)
    # target_link_libraries(pipsipm-unitcommit.exe
    #   ooqpstoch ooqpstochla ooqpmehrotrastoch
    #   ooqpgensparse ooqpbase ooqpsparse ooqpdense
    #   ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})

  if (HAVE_CLANG)
  else()
    add_library(pipsipm-shared SHARED Drivers/pipsipm_C_callbacks.cpp)
    target_link_libraries(pipsipm-shared
      ${WHOLE_ARCHIVE}
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA57_LIBRARY} ${METIS_LIBRARY}
      ${NO_WHOLE_ARCHIVE}
      ${MATH_LIBS})
  endif (HAVE_CLANG)

  # also add MA57 and metis to linear solvers
  set(LINSOLVERS "${LINSOLVERS}${MA57_LIBRARY} ${METIS_LIBRARY}")

endif(HAVE_MA57 AND HAVE_METIS)

#missing linkage mc68 using metis.
#if (HAVE_MA86 AND HAVE_METIS)
#  add_executable(qpgen-sparse-ma86-mehrotra.exe Core/QpGen/QpGenSparseMa86MehrotraDriver.C)
#  target_link_libraries(qpgen-sparse-ma86-mehrotra.exe
#    ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra ooqpdense
#    ${MA86_LIB} ${METIS_LIBRARY} ${MATH_LIBS})
#endif(HAVE_MA86 AND HAVE_METIS)


#build only some executable depending on the linear solvers available
if (HAVE_PARDISO)
  if(HAVE_MA27)
    add_executable(qpgen-sparse-mehrotra.exe Core/QpGen/QpGenSparseMehrotraDriver.C)
    target_link_libraries(qpgen-sparse-mehrotra.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpmehrotra
      ${MA27_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})


  if(HAVE_MA57 AND HAVE_METIS)
    add_executable(qpgen-sparse-ma57-gondzio.exe Core/QpGen/QpGenSparseGondzioDriver.C)
    target_link_libraries(qpgen-sparse-ma57-gondzio.exe
      ooqpgensparse ooqpsparse ooqpbase ooqpgondzio ooqpdense
      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})

    add_executable(pipsipmFromRaw_comm2_schur Drivers/pipsipmFromRaw_comm2_schur.cpp)
    target_link_libraries(pipsipmFromRaw_comm2_schur
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})

    add_executable(pipsipmBatchFromRaw_schur Drivers/pipsipmBatchFromRaw_schur.cpp)
    target_link_libraries(pipsipmBatchFromRaw_schur
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})

    add_executable(pipsipmBatchFromRaw Drivers/pipsipmBatchFromRaw.cpp)
    target_link_libraries(pipsipmBatchFromRaw
      stochInput ${COIN_LIBS}
      ooqpstoch ooqpstochla ooqpmehrotrastoch
      ooqpgensparse ooqpbase ooqpsparse ooqpdense
      ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})

    # add_executable(pipsipmFromRaw_schur32 Drivers/pipsipmFromRaw_schur32.cpp)
    # target_link_libraries(pipsipmFromRaw_schur32
    #   stochInput ${COIN_LIBS}
    #   ooqpstoch ooqpstochla ooqpmehrotrastoch
    #   ooqpgensparse ooqpbase ooqpsparse ooqpdense
    #   ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY32} ${MATH_LIBS})


      add_executable(ooqpFromRaw Drivers/ooqpFromRaw.cpp)
       target_link_libraries(ooqpFromRaw
         stochInput ${COIN_LIBS}
         ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
	 ${MA57_LIBRARY} ${METIS_LIBRARY} ${MA27_LIBRARY} ${MATH_LIBS})

     add_executable(meanSolveAndRecourseEval Drivers/meanSolveAndRecourseEval.cpp)
      target_link_libraries(meanSolveAndRecourseEval
	stochInput ${COIN_LIBS}
	ooqpstoch ooqpstochla ooqpmehrotrastoch
	ooqpgensparse ooqpbase ooqpmehrotra ooqpsparse ooqpdense
	${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY}
	${MATH_LIBS})


#    if(HAVE_AMPL)
#      add_executable(pipsipmPyomo_schur Drivers/pipsipmPyomo_schur.cpp)
#      target_link_libraries(pipsipmPyomo_schur
#	pyomoStochInput stochInput ${COIN_LIBS}
#	ooqpstoch ooqpstochla ooqpmehrotrastoch
#	ooqpgensparse ooqpbase ooqpsparse ooqpdense
#	${AMPL_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${PARDISO_LIBRARY} ${MATH_LIBS})
#    endif(HAVE_AMPL)
  endif(HAVE_MA57 AND HAVE_METIS)
  endif(HAVE_MA27)
endif(HAVE_PARDISO)

#MA27 or MA57 and METIS are prerequisites, this should be fine without additional checks.
#It still works if only one of MA27 or MA57 is available.
add_executable(pipsipmFromRaw Drivers/pipsipmFromRaw.cpp)
target_link_libraries(pipsipmFromRaw
    stochInput ${COIN_LIBS}
    ooqpstoch ooqpstochla ooqpmehrotrastoch
    ooqpgensparse ooqpbase ooqpsparse ooqpdense
    ${MA27_LIBRARY} ${MA57_LIBRARY} ${METIS_LIBRARY} ${MATH_LIBS})
