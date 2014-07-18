cmake_minimum_required (VERSION 2.8.1)
project(linmod)
enable_language(C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


### External dependencies
find_package(R REQUIRED)
include_directories(${R_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)


find_library(LIB_LINMOD NAMES "liblinmod_static.a"
  PATHS ${PROJECT_SOURCE_DIR}/linmod/build NO_DEFAULT_PATH
)


execute_process(COMMAND Rscript "-e" "RNACI:::cppflags(cmake=TRUE)"
                OUTPUT_VARIABLE RNACI_CPPFLAGS)
#execute_process(COMMAND echo ${RNACI_CPPFLAGS})
include_directories(${RNACI_CPPFLAGS})
execute_process(COMMAND Rscript "-e" "RNACI:::ldflags()"
                OUTPUT_VARIABLE RNACI_LDFLAGS)


### Flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfatal-errors -O3")

link_directories(./linmod/build/)



### Source files
set(RLIB_FILES
	glm_fit.c
	glm_utils.c
	lm_fit.c
	utils.c
)



### Make and install
add_library(linmod SHARED ${RLIB_FILES})

target_link_libraries(linmod ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LIB_LINMOD} ${R_LIBRARIES} ${RNACI_LDFLAGS})

set_target_properties(linmod PROPERTIES PREFIX "")
if(NOT WIN32)
  set_target_properties(linmod PROPERTIES SUFFIX ".so")
  install(TARGETS linmod LIBRARY DESTINATION ${CMAKE_SOURCE_DIR})
else(NOT WIN32)
  install(TARGETS linmod RUNTIME DESTINATION ${CMAKE_SOURCE_DIR})
endif(NOT WIN32)

