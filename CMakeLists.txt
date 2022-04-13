cmake_minimum_required(VERSION "3.14")

project(exocross LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" ON)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

add_library(exocross
  accuracy.f90
  input.f90
  phoenix.f90
  spectrum.f90
  timer.f90
  VoigtKampff.f90
)

add_executable(xcross crosssections.f90)
target_link_libraries(xcross exocross)
target_include_directories(xcross PRIVATE ${CMAKE_BINARY_DIR}/modules)
set_target_properties(xcross PROPERTIES PREFIX "")
set_target_properties(xcross PROPERTIES SUFFIX ".exe")

install(TARGETS xcross DESTINATION ${CMAKE_SOURCE_DIR})
