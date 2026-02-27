# @HEADER
# ****************************************************************************
#                Hlangana: Copyright S. Mabuza Enterprises LLC
# 
# Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
# ****************************************************************************
# @HEADER

cmake_minimum_required(VERSION 3.22.0)

project(hlangana)

set(CMAKE_CXX_COMPILER mpicxx)
set(LINKER_LANGUAGE CXX)

set(SOFTWARE "$ENV{SOFTWARE}")
set(PROJECTS "$ENV{PROJECTS}")

find_package(Boost REQUIRED)
if (Boost_FOUND)
  message(STATUS "Found Boost!")
endif()

set(Trilinos_DIR "${SOFTWARE}/trilinos/lib/cmake/Trilinos")
message(STATUS "Looking for Trilinos")
find_package(Trilinos CONFIG REQUIRED)
if( Trilinos_FOUND )
  message(STATUS "Found Trilinos!")
endif()

# Some serial impl in the hydrofem library
set(HYDROFEM_DIR "${SOFTWARE}/hydrofem")
set(HYDROFEM_INCLUDE_DIRS "${HYDROFEM_DIR}/include")
set(HYDROFEM_LIBRARY_DIRS "${HYDROFEM_DIR}/lib")
set(HYDROFEM_LIBRARIES hydrofem)

set(PARMETIS_DIR "${SOFTWARE}/parmetis/")
set(PARMETIS_INCLUDE_DIRS "${SOFTWARE}/parmetis/include")
set(PARMETIS_LIBRARY_DIRS "${SOFTWARE}/parmetis/lib")

set(EIGEN_DIR "${SOFTWARE}/eigen")
set(EIGEN_INCLUDE_DIRS "${EIGEN_DIR}/include/eigen3")
set(EIGEN_LIBRARY_DIRS "${EIGEN_DIR}/lib")

set(CMAKE_CXX_FLAGS "-O3 -lmpi -lparmetis -lmetis -lboost_program_options  -Wno-dev -L${PARMETIS_LIBRARY_DIRS}")

set(CORE "./core")
set(DISC "./disc")
set(MESH "./mesh")
set(LINEAR "./linear")
set(PROJECT_HEADERS "./")

include_directories(${CORE} ${DISC} ${LINEAR} ${PROJECT_HEADERS} ${MESH} ${Trilinos_INCLUDE_DIRS} ${HYDROFEM_INCLUDE_DIRS} ${PARMETIS_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
link_directories(${Trilinos_LIBRARY_DIRS})

file(GLOB SOURCES "${CORE}/*.cpp"
                  "${DISC}/*.cpp"
                  "${MESH}/*.cpp"
                  "${LINEAR}/*.cpp"
                  "${PROJECT_HEADERS}*.cpp")

add_executable(hlangana ${SOURCES})

target_link_libraries(hlangana LINK_PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${HYDROFEM_LIBRARIES} parmetis metis boost_program_options)

install(TARGETS hlangana RUNTIME DESTINATION bin)
