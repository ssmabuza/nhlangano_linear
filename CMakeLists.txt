# ==============================================================================
#                Flujo: Copyright Valiant Scientific
# 
# Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
# ==============================================================================

cmake_minimum_required(VERSION 3.22.0)

project(flujo)

set(CMAKE_CXX_COMPILER mpicxx)
set(LINKER_LANGUAGE CXX)

set(SOFTWARE "$ENV{SOFTWARE}")
set(PROJECTS "$ENV{PROJECTS}")

set(Trilinos_DIR "${SOFTWARE}/trilinos/lib/cmake/Trilinos")
message(STATUS "Looking for Trilinos")
find_package(Trilinos CONFIG REQUIRED)
if( Trilinos_FOUND )
  message(STATUS "Found Trilinos!")
endif()

set(PROJECT_HEADERS "./")

file(GLOB SOURCES "${PROJECT_HEADERS}*.cpp")

add_executable(flujo_mp.exe ${SOURCES})

target_link_libraries(flujo LINK_PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
install(TARGETS flujo RUNTIME DESTINATION bin)
