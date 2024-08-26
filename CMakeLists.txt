cmake_minimum_required(VERSION 3.16)

# Project name and source file language
project(amrex3d LANGUAGES C CXX)

add_executable(amrex3d main.cpp)

if(NOT DEFINED AMReX_ROOT)
  message("-- Download and configure AMReX from Github")

  # Download AMReX from Github
  include(FetchContent)
  set(FETCHCONTENT_QUIET off)

  FetchContent_Declare(
    amrex_code
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git/
    GIT_TAG origin/development
  )

  FetchContent_Populate(amrex_code)

  # CMake will read the files in these directories and configure, build, and install AMReX
  add_subdirectory(${amrex_code_SOURCE_DIR} ${amrex_code_BINARY_DIR})

else ()

  # add AMReX
  message("-- Searching for AMReX install directory at ${AMReX_ROOT}")
  find_package(AMReX PATHS ${AMReX_ROOT}/lib/cmake/AMReX/AMReXConfig.cmake)

  # add mpi
  find_package(MPI REQUIRED)
  target_link_libraries(AMReX::amrex INTERFACE ${MPI_LIBRARIES})
endif()

target_link_libraries(amrex3d PRIVATE AMReX::amrex)


