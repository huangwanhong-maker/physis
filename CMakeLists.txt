# physis: units, geometry, analysis, dynamics, and physical law.
#
# Every module is a target physis::<module>. physis uses calx's and morphe's targets: inside
# The Science Infrastructures, the siblings foundations/calx and foundations/morphe (added
# here when no larger project has added them); built alone elsewhere, they would be fetched
# at pinned versions (not yet written).

cmake_minimum_required(VERSION 3.20)
project(physis VERSION 0.1.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(PHYSIS_TOP_LEVEL ON)
else()
  set(PHYSIS_TOP_LEVEL OFF)
endif()
option(PHYSIS_CHECKS "Build physis's checks" ${PHYSIS_TOP_LEVEL})

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

if(NOT TARGET calx::linalg)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../calx/CMakeLists.txt)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../calx ${CMAKE_CURRENT_BINARY_DIR}/calx)
  else()
    message(FATAL_ERROR "physis needs calx: build it inside The Science Infrastructures (foundations/calx).")
  endif()
endif()
if(NOT TARGET morphe::topology)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../morphe/CMakeLists.txt)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../morphe ${CMAKE_CURRENT_BINARY_DIR}/morphe)
  else()
    message(FATAL_ERROR "physis needs morphe: build it inside The Science Infrastructures (foundations/morphe).")
  endif()
endif()
if(PHYSIS_CHECKS)
  enable_testing()
endif()

# The modules in their order (decision 0005, §2.2).
add_subdirectory(units)
add_subdirectory(analysis)
add_subdirectory(dynamics)
add_subdirectory(quantum)
add_subdirectory(symmetry)
add_subdirectory(geometry)
add_subdirectory(mechanics)
add_subdirectory(spacetime)
add_subdirectory(field_theory)
add_subdirectory(gauge)
add_subdirectory(discretisation)
add_subdirectory(computational_geometry)
add_subdirectory(quantum_gravity)
add_subdirectory(chemistry)
