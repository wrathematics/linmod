cmake_minimum_required(VERSION 2.8.8)
project(linmod)
enable_language(C Fortran)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_Fortran_MODULE_DIRECTORY "include/")

### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfatal-errors -O3 -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wfatal-errors -O3 -mtune=native -fPIC")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### External dependencies
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

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
	src/glm/glm_constants.f90
	src/glm/glm_family_utils.f90
	src/glm/glm_link_utils.f90
	src/glm/glm_loglik_utils.f90
	src/glm/glm_residuals_utils.f90
	src/glm/glm_update_utils.f90
	src/glm/glm.f90
)

set(CINTERFACE_FILES
  src/c_interface/lm_fit.c
)

set(LAPACK_FILES
	src/lapack/lapack.f90
	src/lapack/lapack_omp.f90
)

set(LM_FILES
	src/lm/lm.f90
	src/lm/rdgels.f90
	src/lm/rdgels_utils.f90
	src/lm/rdormqr.f90
	src/lm/rdgeqp3.f90
	src/lm/rdlaqp2.f90
	src/lm/qr_utils.f90
)

set(DISTRIBUTIONS_FILES
	src/distributions/distributions.f90
	src/distributions/cauchy.c
	src/distributions/normal.c
)

set(UTILS_FILES
	src/utils/insertionsorts.f90
	src/utils/quicksort_utils.f90
	src/utils/quicksorts.f90
	src/utils/Rspecial.c
	src/utils/Rspecial_.f90
	src/utils/string_tools.f90
	src/utils/swaps.f90
	src/utils/transposition.f90
	src/utils/linmod_omp.f90
)


set(LINMOD_FILES ${UTILS_FILES} ${DISTRIBUTIONS_FILES} ${LAPACK_FILES} ${LM_FILES} ${GLM_FILES} ${CINTERFACE_FILES})


### Make and install
add_library(linmod_objs OBJECT ${LINMOD_FILES})

add_library(linmod SHARED $<TARGET_OBJECTS:linmod_objs>)
add_library(linmod_static STATIC $<TARGET_OBJECTS:linmod_objs>)


target_link_libraries (linmod ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(linmod PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(linmod PROPERTIES SUFFIX ".so")
  install(TARGETS linmod LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS linmod RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)


#INSTALL(FILES build/libglm.so DESTINATION lib)
