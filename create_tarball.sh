#!/bin/sh

# Cleaning
rm -rf CMakeFiles/ CMakeCache.txt

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
perl -i -pe 's|^add_executable|set\(CMAKE_CXX_FLAGS "\$\{CMAKE_CXX_FLAGS\} -static" CACHE "Flags used by the compiler during tarball build type" STRING\)\nadd_executable|;'   CMakeLists.txt

# Build
cmake .
make

# Strip binary
cp fast fast.unstripped
ls -l fast
strip fast
ls -l fast

# Build tarball
# get fast version
VERSION=`./fast | grep 'FastCodeML V' | head -1 | sed -e 's/^FastCodeML V//'`
rm -rf CMakeFiles/ CMakeCache.txt *.o
DIR=`basename $PWD`
cd ..
mv $DIR FastCodeML-$VERSION
tar cvf FastCodeML-$VERSION.tar --exclude=$(basename $0)  FastCodeML-$VERSION/
gzip -9 FastCodeML-$VERSION.tar
tar tvfz FastCodeML-$VERSION.tar.gz


exit 0

