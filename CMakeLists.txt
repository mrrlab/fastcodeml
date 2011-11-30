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

if(NOT WIN32)
	SET(USE_LAPACK YES CACHE BOOL "Use BLAS/LAPACK")
	SET(USE_DGEMM NO CACHE BOOL "Use DGEMM for UDV computation")
	SET(USE_DSYRK NO CACHE BOOL "Use new method for matrix exponentiation that uses DSYRK")
	SET(USE_MKL_VML NO CACHE BOOL "Use Intel MKL vectorized routines")
endif(NOT WIN32)

SET(USE_NEW_LIKELIHOOD YES CACHE BOOL "Use the new loglikelihood computation method")

SET(USE_OPENMP YES CACHE BOOL "Compile with OpenMP support")
if(USE_OPENMP)
	find_package(OpenMP)
endif(USE_OPENMP)

SET(USE_MPI YES CACHE BOOL "Use MPI for high level parallelization")
if(USE_MPI)
# Find MPI
	#set(MPI_COMPILER "$ENV{MPI_ROOT}/bin/mpicxx")
	#set(MPI_CXX_COMPILER "$ENV{MPI_ROOT}/bin/mpicxx")
	#set(MPI_LIBRARY "$ENV{MPI_ROOT}/lib/libmpich.a")
	#set(MPI_INCLUDE_PATH "$ENV{MPI_ROOT}/include")

	find_package(MPI)
endif(USE_MPI)
if(MPI_FOUND)
	SET(USE_THREAD_MPI YES CACHE BOOL "Use multithread MPI if available")
endif(MPI_FOUND)


# Set options for the various platforms
if($ENV{HOST} MATCHES "matterhorn")

	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
	SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
	SET(CMAKE_SKIP_RPATH)

	add_definitions(-DBOOST_NO_CWCTYPE)
	add_definitions(-xmt2)
	add_definitions(-exceptions)
	add_definitions(-pl ${CMAKE_CURRENT_BINARY_DIR}/fast.pl)
	set(NLOPT_LINK_DIR "${CMAKE_CURRENT_BINARY_DIR}/nlopt-2.2.4" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/nlopt-2.2.4/api" CACHE PATH "NLopt include dir")
	link_directories( ${NLOPT_LINK_DIR} )
	include_directories(${NLOPT_INCLUDE_DIR})
	
	if(USE_NEW_LIKELIHOOD)
		add_definitions(-DNEW_LIKELIHOOD)
	endif(USE_NEW_LIKELIHOOD)

elseif($ENV{HOST} MATCHES "todi")

	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
	SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
	#set(BLAS_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "BLAS lib dir")
	#set(LAPACK_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "LAPACK lib dir")
	#set(MATH_LIB_NAMES acml gfortranbegin gfortran CACHE STRING "Math libraries" )
	#set(MKL_INCLUDE_DIR "/apps/eiger/Intel-FOR-11.1/mkl/include" CACHE PATH "MKL include dir")

	if(USE_LAPACK)
		add_definitions(-DUSE_LAPACK)
	#	link_directories( ${BLAS_LINK_DIR} )
	#	link_directories( ${LAPACK_LINK_DIR} )
		if(USE_DGEMM)
			add_definitions(-DUSE_DGEMM)
		endif(USE_DGEMM)
		if(USE_DSYRK)
			add_definitions(-DUSE_DSYRK)
		endif(USE_DSYRK)
	#	if(USE_MKL_VML)
	#		add_definitions(-DUSE_MKL_VML)
	#	endif(USE_MKL_VML)
	endif(USE_LAPACK)

	set(NLOPT_LINK_DIR "/users/mvalle/apps/lib" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "/users/mvalle/apps/include" CACHE PATH "NLopt include dir")
	link_directories( ${NLOPT_LINK_DIR} )
	include_directories(${NLOPT_INCLUDE_DIR})
	#if(USE_MKL_VML)
	#	include_directories(${MKL_INCLUDE_DIR})
	#endif(USE_MKL_VML)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

	if(MPI_FOUND)
		include_directories(${MPI_INCLUDE_PATH})
		add_definitions(-DUSE_MPI)
		if(USE_THREAD_MPI)
			add_definitions(-DUSE_THREAD_MPI)
		endif(USE_THREAD_MPI)
	endif(MPI_FOUND)
	
	if(USE_NEW_LIKELIHOOD)
		add_definitions(-DNEW_LIKELIHOOD)
	endif(USE_NEW_LIKELIHOOD)

elseif(WIN32)

	set(NLOPT_LINK_DIR "C:/mv/Projects/Selectome/codeml/fast/nlopt/win-build/Debug" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "C:/mv/Projects/Selectome/codeml/fast/nlopt/nlopt-2.2.4/api" CACHE PATH "NLopt include dir")

	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
	add_definitions(/D_SECURE_SCL=0)
	add_definitions(/Oi)
	if(MPI_FOUND)
		include_directories(${MPI_INCLUDE_PATH})
		add_definitions(/DUSE_MPI)
		if(USE_THREAD_MPI)
			add_definitions(/DUSE_THREAD_MPI)
		endif(USE_THREAD_MPI)
	endif(MPI_FOUND)
	
	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(OPENMP_FOUND)

	link_directories( ${NLOPT_LINK_DIR} )
	include_directories(${NLOPT_INCLUDE_DIR})
	
	if(USE_NEW_LIKELIHOOD)
		add_definitions(/DNEW_LIKELIHOOD)
	endif(USE_NEW_LIKELIHOOD)

else()

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

	else(CMAKE_COMPILER_IS_GNUCXX)
		SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
		SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
		SET(CMAKE_SKIP_RPATH)
	endif(CMAKE_COMPILER_IS_GNUCXX)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

	if(MPI_FOUND)
		include_directories(${MPI_INCLUDE_PATH})
		add_definitions(-DUSE_MPI)
		if(USE_THREAD_MPI)
			add_definitions(-DUSE_THREAD_MPI)
		endif(USE_THREAD_MPI)
	endif(MPI_FOUND)

	set(BLAS_LINK_DIR "/apps/eiger/Intel-CPP-11.1/mkl/lib/em64t" CACHE PATH "BLAS lib dir")
	set(LAPACK_LINK_DIR "/apps/eiger/Intel-CPP-11.1/mkl/lib/em64t" CACHE PATH "LAPACK lib dir")
	#set(MATH_LIB_NAMES blas lapack gfortranbegin gfortran CACHE STRING "Math libraries" )
	set(MATH_LIB_NAMES mkl_intel_lp64 mkl_sequential mkl_lapack mkl_core m CACHE STRING "Math libraries" )
	set(MKL_INCLUDE_DIR "/apps/eiger/Intel-CPP-11.1/mkl/include" CACHE PATH "MKL include dir")

	if(USE_LAPACK)
		add_definitions(-DUSE_LAPACK)
		link_directories( ${BLAS_LINK_DIR} )
		link_directories( ${LAPACK_LINK_DIR} )
		if(USE_DGEMM)
			add_definitions(-DUSE_DGEMM)
		endif(USE_DGEMM)
		if(USE_DSYRK)
			add_definitions(-DUSE_DSYRK)
		endif(USE_DSYRK)
		if(USE_MKL_VML)
			add_definitions(-DUSE_MKL_VML)
		endif(USE_MKL_VML)
	endif(USE_LAPACK)

	set(NLOPT_LINK_DIR "/users/mvalle/apps/lib" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "/users/mvalle/apps/include" CACHE PATH "NLopt include dir")
	link_directories( ${NLOPT_LINK_DIR} )
	include_directories(${NLOPT_INCLUDE_DIR})
	if(USE_MKL_VML)
		include_directories(${MKL_INCLUDE_DIR})
	endif(USE_MKL_VML)
	
	if(USE_NEW_LIKELIHOOD)
		add_definitions(-DNEW_LIKELIHOOD)
	endif(USE_NEW_LIKELIHOOD)

endif()


# Executable
add_executable(fast ${SRCS})


# Link libraries
if(WIN32)
	target_link_libraries(fast nlopt)
else(WIN32)
	if(USE_LAPACK)
		target_link_libraries(fast ${MATH_LIB_NAMES} nlopt m)
	else(USE_LAPACK)
		target_link_libraries(fast nlopt m)
	endif(USE_LAPACK)
endif(WIN32)

if(MPI_FOUND)
	target_link_libraries(fast ${MPI_LIBRARY})
endif(MPI_FOUND)

# Make documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
   set(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")

  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG})
endif(DOXYGEN_FOUND)

