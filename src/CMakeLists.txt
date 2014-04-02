cmake_minimum_required (VERSION 2.8.1)
project(linmod)
enable_language(C Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})


execute_process(COMMAND Rscript "-e" "SEXPtools:::cppflags(cmake=TRUE)"
                OUTPUT_VARIABLE SEXPTOOLS_CPPFLAGS)
#execute_process(COMMAND echo ${SEXPTOOLS_CPPFLAGS})
include_directories(${SEXPTOOLS_CPPFLAGS})
execute_process(COMMAND Rscript "-e" "SEXPtools:::ldflags()"
                OUTPUT_VARIABLE SEXPTOOLS_LDFLAGS)


find_package(OpenMP)

if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  
  # OpenMP flags for Fortran --- thanks internet man who I stole this from!
  # http://sethmmorton.github.io/How%20To/2013/04/28/using-cmake-to-find-openmp-flags-for-fortran-compilers/
  find_package(OpenMP_Fortran)
  set_target_properties(${EXECUTABLE_NAME} PROPERTIES
                      COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                      LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
  
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfatal-errors -O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Wimplicit-interface")



### Source files
set(GLM_FILES
	linmod/src/glm/glm_check.f90
	linmod/src/glm/glm_link_utils.f90
	linmod/src/glm/glm_loglik_utils.f90
	linmod/src/glm/glm_mu_var.f90
	linmod/src/glm/glm_residuals_utils.f90
	linmod/src/glm/glm_update_utils.f90
	linmod/src/glm/glm.f90
)

set(CINTERFACE_FILES
  linmod/src/c_interface/lm_fit.c
)

set(LAPACK_FILES
	linmod/src/lapack/dgeadd_omp.f90
	linmod/src/lapack/dlacpy.f90
	linmod/src/lapack/lapack.f90
	linmod/src/lapack/lapack_omp.f90
)

set(LM_FILES
	linmod/src/lm/lmfit.f90
	linmod/src/lm/rdgels.f90
	linmod/src/lm/rdgeqpf.f90
	linmod/src/lm/rdormqr.f90
	linmod/src/lm/rdgeqp3.f90
	linmod/src/lm/rdlaqp2.f90
)

set(UTILS_FILES
	linmod/src/utils/string_tools.f90
)

set(RLIB_FILES
	glm_fit.c
	lm_fit.c
)



### Make and install
add_library(linmod SHARED ${UTILS_FILES} ${LAPACK_FILES} ${LM_FILES} ${GLM_FILES} ${CINTERFACE_FILES} ${RLIB_FILES})

target_link_libraries (linmod ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${R_LIBRARIES} ${SEXPTOOLS_LDFLAGS})

set_target_properties(linmod PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(linmod PROPERTIES SUFFIX ".so")
  install(TARGETS linmod LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS linmod RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)

