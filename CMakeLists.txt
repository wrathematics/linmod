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
  src/glm.f90
  src/glm_link_utils.f90
  src/glm_loglik_utils.f90
  src/glm_mu_var.f90
)



add_library(glm SHARED ${GLM_FILES})
add_library(glm_static STATIC ${GLM_FILES})


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
