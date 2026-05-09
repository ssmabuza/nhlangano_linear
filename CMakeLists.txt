# ==============================================================================
#                Flujo: Copyright Valiant Scientific
# 
# Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
# ==============================================================================

cmake_minimum_required(VERSION 3.22.0)

project(flujo)

message(STATUS "Looking for Trilinos")
find_package(Trilinos CONFIG REQUIRED)
if( Trilinos_FOUND )
  message(STATUS "Found Trilinos!")
endif()

# Use the same compiler as Trilinos
set(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)
set(LINKER_LANGUAGE CXX)
set(CMAKE_CXX_STANDARD 20)

add_executable(flujo
  Flujo_Main.cpp
  src/driver/Flujo_Driver.cpp
)

target_include_directories(flujo PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/driver"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/closures"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/eqn_sets"
)

target_link_libraries(flujo LINK_PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
install(TARGETS flujo RUNTIME DESTINATION bin)
