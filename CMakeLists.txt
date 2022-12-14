cmake_minimum_required(VERSION 2.8.0)
project(FastCodeML)

# Source code
set(SRCS
	fast.cpp
	CmdLine.cpp
	Genes.cpp
	Phylip.cpp
	PhyloTree.cpp
	Newick.cpp
	TreeNode.cpp
	BayesTest.cpp
	FillMatrix.cpp
	Forest.cpp
	TransitionMatrix.cpp
	BranchSiteModel.cpp
	ProbabilityMatrixSet.cpp
	FatVectorTransform.cpp
	CodonFrequencies.cpp
	AlignedAllocator.cpp
	HighLevelCoordinator.cpp
	CodeMLoptimizer.cpp
	ForestExport.cpp
	ParseParameters.cpp
	VerbosityLevels.cpp
	DAGScheduler.cpp
	TreeAndSetsDependencies.cpp
	WriteResults.cpp
)


# Search for Boost version 1.55 minimum
set(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.55")
find_package(Boost 1.42)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message(SEND_ERROR "Boost::spirit required")
endif(Boost_FOUND)


# Get the *basic* configuration switches
OPTION(USE_LAPACK         		"Use BLAS/LAPACK" OFF)
OPTION(USE_OPENMP         		"Compile with OpenMP support" ON)


# Get the *advanced* configuration switches (change default values with care, some options are not fully functional)
OPTION(USE_MKL_VML        		"Use Intel MKL vectorized routines" OFF) 
OPTION(USE_MPI            		"Use MPI for high level parallelization" OFF)
if(NOT WIN32)
	OPTION(BUILD_NOT_SHARED   	"Build FastCodeML not shared" OFF) 
endif(NOT WIN32)
OPTION(BUILD_SEARCH_MPI   		"Search for MPI installation?" OFF) 
OPTION(USE_ORIGINAL_PROPORTIONS		"Use the original CodeML proportion definition" ON)
SET(USE_LIKELIHOOD_METHOD		"Original" CACHE STRING "Select the type of likelihood computation method: Original, NonRecursive, FatVector, DAG")
SET_PROPERTY(CACHE USE_LIKELIHOOD_METHOD PROPERTY STRINGS Original NonRecursive FatVector DAG)
OPTION(USE_IDENTITY_MATRIX		"Force identity matrix when time is zero" OFF)
OPTION(USE_CPV_SCALING			"Scale conditional probability vectors to avoid under/overflow" ON)


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
	#add_definitions(/arch:AVX)
else(WIN32)
	if(BUILD_NOT_SHARED)
		SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
		SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
		SET(CMAKE_SKIP_RPATH ON)
	endif(BUILD_NOT_SHARED)
endif(WIN32)

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

if(USE_IDENTITY_MATRIX)
	add_definitions(-DFORCE_IDENTITY_MATRIX)
endif(USE_IDENTITY_MATRIX)

# Select on method to compute likelihood
if(USE_LIKELIHOOD_METHOD MATCHES "NonRecursive")
	add_definitions(-DNON_RECURSIVE_VISIT)
	set(USE_CPV_SCALING OFF)
endif(USE_LIKELIHOOD_METHOD MATCHES "NonRecursive")

if(USE_LIKELIHOOD_METHOD MATCHES "FatVector")
	add_definitions(-DNEW_LIKELIHOOD)
	set(USE_CPV_SCALING OFF)
endif(USE_LIKELIHOOD_METHOD MATCHES "FatVector")

if(USE_LIKELIHOOD_METHOD MATCHES "DAG")
	add_definitions(-DUSE_DAG)
	set(USE_CPV_SCALING OFF)
endif(USE_LIKELIHOOD_METHOD MATCHES "DAG")

# CPV scaling is valid only for the Original likelihood method
if(USE_CPV_SCALING)
	add_definitions(-DUSE_CPV_SCALING)
endif(USE_CPV_SCALING)

# Set paths
if(USE_MPI)
	if(MPI_C_INCLUDE_PATH)
		include_directories(${MPI_C_INCLUDE_PATH})
	endif(MPI_C_INCLUDE_PATH)
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

###########################################################
if(NOT WIN32)

#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_CXX_FLAGS_DEBUG "-g -error -pg -O -Wall -std=c++98 -Wextra -Wno-unused-parameter -Wunused" CACHE "Debug mode options" STRING)
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb -g -pg -O0 -Wall -std=c++98 -Wextra -Wno-unused-parameter -Wunused")
endif(NOT WIN32)

###########################################################

# Executable
add_executable(fast ${SRCS})

# Link libraries
target_link_libraries(fast nlopt)
if(MATH_LIB_NAMES AND USE_LAPACK)
	target_link_libraries(fast ${MATH_LIB_NAMES})
endif(MATH_LIB_NAMES AND USE_LAPACK)

target_link_libraries(fast ${EXTRA_LIBS})

if(MPI_LIBRARY)
	target_link_libraries(fast ${MPI_LIBRARY})
endif(MPI_LIBRARY)