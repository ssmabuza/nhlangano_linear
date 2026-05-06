#!/bin/bash

TRILINOS_HOME="${HOME}/Software/trilinos/Trilinos"

INSTALL_DIR="${HOME}/Software/trilinos"

SOFTWARE="${HOME}/Software"

MPI_DIR="${MPI_DIR}/bin/mpi"

#      -D CMAKE_CXX_FLAGS_DEBUG_OVERRIDE:STRING="-g -O3 -lgfortran -ansi -pedantic -Wall -Wno-long-long -Wno-unused-local-typedefs -Wno-strict-aliasing -DBOOST_NO_HASH -Wno-virtual-move-assign"\
#       -D CMAKE_CXX_FLAGS:STRING="-g -O3 -fPIC -funroll-loops -ansi -pedantic -ftrapv -Wall -Wno-long-long -Wno-strict-aliasing -DBOOST_NO_HASH -DSKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS" \
#       -D Trilinos_EXTRA_LINK_FLAGS:STRING="-L/usr/lib64 /usr/lib/x86_64-linux-gnu/libcurl.so.4 /usr/lib/x86_64-linux-gnu/libz.so -lhdf5_hl -lhdf5 -ldl" \


rm -rf CMakeCache.txt CMakeFiles 
\
cmake \
      -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
      -D Trilinos_ENABLE_INSTALL_CMAKE_CONFIG_FILES:BOOL=ON \
      -D CMAKE_CXX_STANDARD=20 \
      -D CMAKE_CXX_STANDARD_REQUIRED=ON \
      -D CMAKE_CXX_EXTENSIONS=OFF \
      -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
      -D Trilinos_ENABLE_TESTS:BOOL=OFF \
      -D Kokkos_ENABLE_TESTS:BOOL=OFF \
      -D BUILD_SHARED_LIBS:BOOL=ON \
      -D Trilinos_ENABLE_DEBUG=OFF \
      -D Trilinos_ENABLE_OpenMP=OFF \
      -D Kokkos_ENABLE_Pthread:BOOL=OFF \
      -D Tpetra_INST_INT_INT:BOOL=ON \
      -D CMAKE_BUILD_TYPE:STRING=RELEASE \
      -D Phalanx_KOKKOS_DEVICE_TYPE:STRING="SERIAL" \
      -D Trilinos_ENABLE_Fortran:BOOL=ON \
      -D CMAKE_C_FLAGS_DEBUG_OVERRIDE:STRING="-g -O3" \
      -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
      -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \
      -D Trilinos_ENABLE_Teko:BOOL=ON \
      -D Teko_ENABLE_TESTS:BOOL=OFF \
      -D Teko_ENABLE_EXAMPLES:BOOL=OFF \
      -D Trilinos_ENABLE_MueLu:BOOL=ON \
      -D MueLu_ENABLE_Experimental:BOOL=ON \
      -D Trilinos_ENABLE_Belos:BOOL=ON \
      -D Trilinos_ENABLE_Panzer:BOOL=ON \
      -D Trilinos_ENABLE_Shards:BOOL=ON \
      -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
      -D Trilinos_ENABLE_Zoltan:BOOL=ON \
      -D Trilinos_ENABLE_Amesos2:BOOL=ON \
      -D Trilinos_ENABLE_SEACAS:BOOL=ON \
      -D Trilinos_ENABLE_SEACASIoss:BOOL=ON \
      -D Trilinos_ENABLE_STK:BOOL=ON \
      -D Trilinos_ENABLE_STKMesh:BOOL=ON \
      -D Trilinos_ENABLE_STKUtil:BOOL=ON \
      -D Trilinos_ENABLE_STKTopology:BOOL=ON \
      -D Trilinos_ENABLE_STKIO:BOOL=ON \
      -D Trilinos_ENABLE_STKTransfer:BOOL=ON \
      -D Trilinos_ENABLE_STKSearchUtil:BOOL=ON \
      -D Trilinos_ENABLE_STKSearch:BOOL=ON \
      -D Trilinos_ENABLE_STKUnit_test_utils:BOOL=OFF \
      -D Trilinos_ENABLE_Stokhos:BOOL=OFF \
      -D Trilinos_ENABLE_Xpetra:BOOL=ON \
      -D Xpetra_ENABLE_Experimental:BOOL=ON \
      -D Panzer_ENABLE_TESTS:BOOL=ON \
      -D Panzer_ENABLE_EXAMPLES:BOOL=OFF \
      -D Panzer_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
      -D Trilinos_ENABLE_PanzerAdaptersSTK:BOOL=ON \
      -D SEACASExodus_ENABLE_MPI:BOOL=ON \
      -D Trilinos_ENABLE_STKTools:BOOL=ON \
      -D Intrepid2_ENABLE_DEBUG_INF_CHECK:BOOL=OFF \
      -D CMAKE_CXX_COMPILER:FILEPATH="${MPI_DIR}/bin/mpicxx" \
      -D CMAKE_C_COMPILER:FILEPATH="${MPI_DIR}/bin/mpicc" \
      -D CMAKE_Fortran_COMPILER:FILEPATH="${MPI_DIR}/bin/mpif90" \
      -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
      -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
      -D CMAKE_INSTALL_PREFIX:PATH="${INSTALL_DIR}" \
      -D Trilinos_VERBOSE_CONFIGURE:BOOL=OFF \
      -D Trilinos_ENABLE_STRONG_CXX_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_STRONG_C_COMPILE_WARNINGS=OFF \
      -D Trilinos_ENABLE_SHADOW_WARNINGS=OFF \
      -D TPL_ENABLE_MPI:BOOL=ON \
      -D MPI_BASE_DIR:PATH="${MPI_DIR}" \
      -D TPL_ENABLE_Boost:BOOL=OFF \
      -D TPL_ENABLE_BoostLib:BOOL=OFF \
      -D TPL_ENABLE_Netcdf:BOOL=ON \
      -D TPL_Netcdf_INCLUDE_DIRS:PATH="${SOFTWARE}/trilinos/include" \
      -D TPL_Netcdf_LIBRARIES:FILEPATH="${SOFTWARE}/trilinos/lib/libnetcdf.so" \
      -D TPL_ENABLE_HDF5:BOOL=ON \
      -D HDF5_INCLUDE_DIRS:PATH="${SOFTWARE}/trilinos/include" \
      -D HDF5_LIBRARY_DIRS:PATH="${SOFTWARE}/trilinos/lib" \
      -D TPL_ENABLE_LAPACK:BOOL=ON \
      -D TPL_LAPACK_LIBRARIES:FILEPATH="/usr/lib/x86_64-linux-gnu/liblapack.so" \
      -D TPL_ENABLE_BLAS:BOOL=ON \
      -D TPL_BLAS_LIBRARIES:FILEPATH="/usr/lib/x86_64-linux-gnu/libblas.so" \
      -D TPL_ENABLE_gtest:BOOL=OFF \
      -D CMAKE_SKIP_RULE_DEPENDENCY=ON \
      -D TPL_ENABLE_Matio:BOOL=OFF \
      -D TPL_ENABLE_X11:BOOL=OFF \
      -D TPL_ENABLE_SuperLU:BOOL=OFF \
${TRILINOS_HOME}
