#!/bin/sh

# Cleaning
rm -rf CMakeFiles/ CMakeCache.txt Makefile cmake_install.cmake fast.unstripped fast
dos2unix CMakeLists.txt

# Set variables used at compile time, assume Vital-IT installation paths
export PATH=/software/bin${PATH:+:$PATH}
export BLAS_LIB_DIR="/software/Utility/OpenBLAS/0.2.5/lib"
export LAPACK_LIB_DIR="/software/Utility/OpenBLAS/0.2.5/include"
export NLOPT_LIB_DIR="/software/Utility/nlopt/2.3/lib"
export NLOPT_INCLUDE_DIR="/software/Utility/nlopt/2.3/include"
export MATH_LIB_NAMES="openblas;lapack;gfortranbegin;gfortran"
export CXX=g++

# Set right options in CMakeLists.txt file
perl -i -pe 's|^(OPTION\(USE_LAPACK .*) OFF\)$|$1 ON\)|;'                   CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_MKL_VML .*) ON\)$|$1 OFF\)|;'                  CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_OPENMP .*) OFF\)$|$1 ON\)|;'                   CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_MPI .*) ON\)$|$1 OFF\)|;'                      CMakeLists.txt
perl -i -pe 's|^(OPTION\(BUILD_SEARCH_MPI .*) ON\)$|$1 OFF\)|;'             CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_ORIGINAL_PROPORTIONS .*) ON\)$|$1 OFF\)|;'     CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_IDENTITY_MATRIX .*) ON\)$|$1 OFF\)|;'          CMakeLists.txt
perl -i -pe 's|^(OPTION\(USE_CPV_SCALING .*) OFF\)$|$1 ON\)|;'              CMakeLists.txt

# Set right compilation options (remove debugging options)
perl -i -pe 's|^set\(CMAKE_CXX_FLAGS_DEBUG|#set\(CMAKE_CXX_FLAGS_DEBUG|;'   CMakeLists.txt
perl -i -pe 's|^add_executable|SET\(CMAKE_BUILD_TYPE RELEASE\)\nset\(CMAKE_EXE_LINKER_FLAGS_RELEASE "-static" CACHE "Release mode linker options" STRING FORCE\)\nadd_executable|;'   CMakeLists.txt

# Build
cmake .
make

# Strip binary
cp fast fast.unstripped
ls -l fast
strip fast
ls -l fast

# Build tarball
VERSION=`./fast | grep 'FastCodeML V' | head -1 | sed -e 's/^FastCodeML V//'`
rm -rf CMakeFiles/ CMakeCache.txt Makefile cmake_install.cmake *.o CMakeLists.txt
svn up
mkdir FastCodeML-$VERSION
cp -r * FastCodeML-$VERSION/
rm -Rf FastCodeML-$VERSION/FastCodeML-$VERSION/
tar cvf FastCodeML-$VERSION.tar --exclude=$(basename $0) --exclude=.svn --exclude=fast.unstripped --exclude=TODO  FastCodeML-$VERSION/
gzip -9 FastCodeML-$VERSION.tar
tar tvfz FastCodeML-$VERSION.tar.gz
rm -Rf FastCodeML-$VERSION/


exit 0

