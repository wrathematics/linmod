cmake_minimum_required (VERSION 2.6)
project(glm)
enable_language(C Fortran)

### Flags
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -w")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -w")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

#find_package(Check)
#find_package(Cobra)

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


### Source files
set(GLM_FILES
	src/glm/glm_check.f90
	src/glm/glm_link_utils.f90
	src/glm/glm_loglik_utils.f90
	src/glm/glm_mu_var.f90
	src/glm/glm_residuals_utils.f90
	src/glm/glm_update_utils.f90
	src/glm/glm.f90
)

set(CINTERFACE_FILES
  src/c_interface/lm_fit.c
)

set(LAPACK_FILES
	src/lapack/dgeadd_omp.f90
	src/lapack/dlacpy.f90
	src/lapack/lapack.f90
	src/lapack/lapack_omp.f90
)

set(LM_FILES
	src/lm/lmfit.f90
	src/lm/rdgels.f90
	src/lm/rdgeqpf.f90
	src/lm/rdormqr.f90
)

set(UTILS_FILES
	src/utils/string_tools.f90
)


add_library(glm SHARED ${UTILS_FILES} ${LAPACK_FILES} ${LM_FILES} ${GLM_FILES} ${CINTERFACE_FILES})
add_library(glm_static STATIC ${UTILS_FILES} ${LAPACK_FILES} ${LM_FILES} ${GLM_FILES} ${CINTERFACE_FILES})


target_link_libraries (glm ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

#if (CHECK_FOUND OR COBRA_FOUND)
#	enable_testing()

#	if(CHECK_FOUND)
#		add_subdirectory(test)
#	endif()

#	if(COBRA_FOUND)
#		add_subdirectory(bench)
#	endif()
#endif()

INSTALL(FILES build/libglm.so DESTINATION lib)