cmake_minimum_required(VERSION 2.6)
project(FastCodeML)

# Source code
set(SRCS
	fast.cpp
	CmdLine.cpp
	Genes.cpp
	PhyloTree.cpp
	TreeNode.cpp
	BayesTest.cpp
	FillQ.cpp
	Forest.cpp
	TransitionMatrix.cpp
	BranchSiteModel.cpp
	TransitionMatrixSet.cpp
	FatVectorTransform.cpp
	CodonFrequencies.cpp
	AlignedAllocator.cpp
	HighLevelCoordinator.cpp
)


# Search for Boost version 1.42 minimum
set(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0")
find_package(Boost 1.42)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message(SEND_ERROR "Boost::spirit required")
endif(Boost_FOUND)


# Get the configuration switches
OPTION(USE_LAPACK         "Use BLAS/LAPACK" ON)
OPTION(USE_MKL_VML        "Use Intel MKL vectorized routines" OFF)
OPTION(USE_NEW_LIKELIHOOD "Use the new loglikelihood computation method" OFF)
OPTION(USE_OPENMP         "Compile with OpenMP support" ON)
OPTION(USE_MPI            "Use MPI for high level parallelization" ON)
if(NOT WIN32)
OPTION(BUILD_NOT_SHARED   "Build FastCodeML not shared" OFF)
endif(NOT WIN32)
OPTION(BUILD_SEARCH_MPI   "Search for MPI installation?" OFF)
OPTION(USE_ORIGINAL_PROPORTIONS "Use the original CodeML proportion definition" OFF)


# Search for OpenMP support
if(USE_OPENMP)
	find_package(OpenMP)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		if(NOT WIN32)
			set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
		endif(NOT WIN32)
	endif(OPENMP_FOUND)

endif(USE_OPENMP)


# Search for MPI
if(USE_MPI)
	if(BUILD_SEARCH_MPI)
		find_package(MPI)
	else(BUILD_SEARCH_MPI)
		set(MPI_FOUND YES)
	endif(BUILD_SEARCH_MPI)
endif(USE_MPI)


# Load the needed paths
set(LINK_DIR_BLAS     $ENV{BLAS_LIB_DIR}      CACHE PATH "BLAS lib dir")
set(LINK_DIR_LAPACK   $ENV{LAPACK_LIB_DIR}    CACHE PATH "LAPACK lib dir")
set(INCLUDE_DIR_MKL   $ENV{MKL_INCLUDE_DIR}   CACHE PATH "MKL include dir")
set(INCLUDE_DIR_NLOPT $ENV{NLOPT_INCLUDE_DIR} CACHE PATH "NLopt include dir")
set(LINK_DIR_NLOPT    $ENV{NLOPT_LIB_DIR}     CACHE PATH "NLopt lib dir")
set(MATH_LIB_NAMES    $ENV{MATH_LIB_NAMES}    CACHE STRING "Math libraries (Separated by ';')") 


# Set compiler switches
if(WIN32)
	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
	add_definitions(/D_SECURE_SCL=0)
	add_definitions(/Oi)
else(WIN32)
	if(BUILD_NOT_SHARED)
		SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
		SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
		SET(CMAKE_SKIP_RPATH)
	endif(BUILD_NOT_SHARED)
endif(WIN32)

if(USE_NEW_LIKELIHOOD)
	add_definitions(-DNEW_LIKELIHOOD)
endif(USE_NEW_LIKELIHOOD)

if(USE_MKL_VML)
	add_definitions(-DUSE_MKL_VML)
endif(USE_MKL_VML)

if(USE_LAPACK)
	add_definitions(-DUSE_LAPACK)
endif(USE_LAPACK)

if(MPI_FOUND)
	add_definitions(-DUSE_MPI)
endif(MPI_FOUND)

if(USE_ORIGINAL_PROPORTIONS)
	add_definitions(-DUSE_ORIGINAL_PROPORTIONS)
endif(USE_ORIGINAL_PROPORTIONS)


# Set paths
if(USE_MPI)
	if(MPI_INCLUDE_PATH)
		include_directories(${MPI_INCLUDE_PATH})
	endif(MPI_INCLUDE_PATH)
	if(MPI_LIBRARY)
		link_directories(${MPI_LIBRARY})
	endif(MPI_LIBRARY)
endif(USE_MPI)
if(USE_MKL_VML)
	include_directories(${INCLUDE_DIR_MKL})
endif(USE_MKL_VML)
if(USE_LAPACK)
	link_directories(${LINK_DIR_BLAS})
	link_directories(${LINK_DIR_LAPACK})
endif(USE_LAPACK)
link_directories(${LINK_DIR_NLOPT})
include_directories(${INCLUDE_DIR_NLOPT})


# Set options for the various platforms
if($ENV{HOST} MATCHES "matterhorn")

	add_definitions(-DBOOST_NO_CWCTYPE)
	add_definitions(-xmt2)
	add_definitions(-exceptions)
	add_definitions(-pl ${CMAKE_CURRENT_BINARY_DIR}/fast.pl)

endif($ENV{HOST} MATCHES "matterhorn")

if($ENV{HOST} MATCHES "julier")

	set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-mkl" CACHE "Release mode linker options" STRING FORCE)
	set(CMAKE_CXX_FLAGS_RELEASE "-fast -xSSE4.2 -mkl -restrict -fstrict-aliasing -Wstrict-aliasing -fomit-frame-pointer -unroll -finline-functions" CACHE "Release mode options" STRING FORCE)

endif($ENV{HOST} MATCHES "julier")

if($ENV{HOST} MATCHES "eiger")

	#Set g++ specific compiler settings
	if(CMAKE_COMPILER_IS_GNUCXX)

		#slow, but suitable for debugging
		if(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES "-Wextra")
			set(CMAKE_CXX_FLAGS_DEBUG "-g -O -Wall -std=c++98 -Wextra -Wno-unused-parameter" CACHE "Debug mode options" STRING FORCE)
			#  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wno-unused-result -Wno-write-strings -g -pg -fcheck-data-deps")
			#  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -pg -fcheck-data-deps" CACHE STRING "Debug mode options")
		endif(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES "-Wextra")

		# fast; -fomit-frame-pointer and -finline-functions not needed, automatically set by -O3; -O4 or higher not supported by gcc;
		# -funroll-loops might be useful, but can sometimes make the code slower (see man gcc). We suppress warnings regarding unused results.
		# -O3 -Wall -std=c++98 -Wextra -pedantic -ffast-math -mtune=native -minline-stringops-dynamically -Wno-unused-result -funroll-loops -funsafe-loop-optimizations -Wunsafe-loop-optimizations"

		if(NOT ${CMAKE_CXX_FLAGS_RELEASE} MATCHES "-minline-stringops-dynamically")
			set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fstrict-aliasing -Wstrict-aliasing -std=c++98 -Wextra -ffast-math -msse2 -mtune=native -minline-stringops-dynamically -funroll-loops -fbuiltin" CACHE "Release mode options" STRING FORCE)
		endif(NOT ${CMAKE_CXX_FLAGS_RELEASE} MATCHES "-minline-stringops-dynamically")

		#if(NOT ${CMAKE_EXE_LINKER_FLAGS_RELEASE} MATCHES "Wl")
			set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-Wl,-O1" CACHE "Release mode linker options" STRING FORCE)
		#endif(NOT ${CMAKE_EXE_LINKER_FLAGS_RELEASE} MATCHES "Wl")
	endif(CMAKE_COMPILER_IS_GNUCXX)

endif($ENV{HOST} MATCHES "eiger")


# Executable
add_executable(fast ${SRCS})


# Link libraries
target_link_libraries(fast nlopt)
if(MATH_LIB_NAMES)
	target_link_libraries(fast ${MATH_LIB_NAMES})
endif(MATH_LIB_NAMES)

if(MPI_LIBRARY)
	target_link_libraries(fast ${MPI_LIBRARY})
endif(MPI_LIBRARY)


# Make documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
   set(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")

  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG})
endif(DOXYGEN_FOUND)

