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
)

# Search for Boost version 1.43
set(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0")
find_package(Boost 1.42)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
else(Boost_FOUND)
	message(SEND_ERROR "Boost::spirit required")
endif(Boost_FOUND)

# Select build options
SET(USE_OPTIMIZER YES CACHE BOOL "Use the optimizer")
if(NOT WIN32)
	SET(USE_LAPACK YES CACHE BOOL "Use BLAS/LAPACK")
endif(NOT WIN32)

SET(USE_OPENMP YES CACHE BOOL "Compile with OpenMP")
if(USE_OPENMP)
	find_package(OpenMP)
endif(USE_OPENMP)

# Set BLAS/LAPACK directories and libraries
if($ENV{HOST} MATCHES "rosa")

	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
	SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
	set(BLAS_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "BLAS lib dir")
	set(LAPACK_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "LAPACK lib dir")
	set(MATH_LIB_NAMES acml gfortranbegin gfortran CACHE STRING "Math libraries" )

	#add_definitions(-DFORCE_ARGUMENTS)

	if(USE_LAPACK)
		add_definitions(-DUSE_LAPACK)
		link_directories( ${BLAS_LINK_DIR} )
		link_directories( ${LAPACK_LINK_DIR} )
	endif(USE_LAPACK)

	set(NLOPT_LINK_DIR "/users/mvalle/apps/lib" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "/users/mvalle/apps/include" CACHE PATH "NLopt include dir")
	if(USE_OPTIMIZER)
		link_directories( ${NLOPT_LINK_DIR} )
		include_directories(${NLOPT_INCLUDE_DIR})
		add_definitions(-DUSE_OPTIMIZER)
	endif(USE_OPTIMIZER)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

elseif($ENV{HOST} MATCHES "palu")

	SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
	SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
	set(BLAS_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "BLAS lib dir")
	set(LAPACK_LINK_DIR "/opt/acml/4.4.0/pgi64/lib" CACHE PATH "LAPACK lib dir")
	set(MATH_LIB_NAMES acml gfortranbegin gfortran CACHE STRING "Math libraries" )

	#add_definitions(-DFORCE_ARGUMENTS)

	if(USE_LAPACK)
		add_definitions(-DUSE_LAPACK)
		link_directories( ${BLAS_LINK_DIR} )
		link_directories( ${LAPACK_LINK_DIR} )
	endif(USE_LAPACK)

	set(NLOPT_LINK_DIR "/users/mvalle/apps/lib" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "/users/mvalle/apps/include" CACHE PATH "NLopt include dir")
	if(USE_OPTIMIZER)
		link_directories( ${NLOPT_LINK_DIR} )
		include_directories(${NLOPT_INCLUDE_DIR})
		add_definitions(-DUSE_OPTIMIZER)
	endif(USE_OPTIMIZER)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

elseif(WIN32)

	set(NLOPT_LINK_DIR "C:/mv/Projects/Selectome/codeml/fast/nlopt/Release" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "C:/mv/Projects/Selectome/codeml/fast/nlopt/nlopt-2.2.1-pc/api" CACHE PATH "NLopt include dir")

	add_definitions(/D_CRT_SECURE_NO_WARNINGS)
	add_definitions(/D_SECURE_SCL=0)
	add_definitions(/Oi)
	#add_definitions(/DFORCE_ARGUMENTS)
	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

	if(USE_OPTIMIZER)
		link_directories( ${NLOPT_LINK_DIR} )
		include_directories(${NLOPT_INCLUDE_DIR})
		add_definitions(/DUSE_OPTIMIZER)
	endif(USE_OPTIMIZER)

else()

	#Set g++ specific compiler settings
	if(CMAKE_COMPILER_IS_GNUCXX)

		#slow, but suitable for debugging
		if(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES "-Wextra")
			set(CMAKE_CXX_FLAGS_DEBUG "-g -O -Wall -std=c++98 -Wextra -Wno-unused-parameter -pedantic -Wshadow" CACHE "Debug mode options" STRING FORCE)
			#  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -Wno-unused-result -Wno-write-strings -g -pg -fcheck-data-deps")
			#  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -pg -fcheck-data-deps" CACHE STRING "Debug mode options")
		endif(NOT ${CMAKE_CXX_FLAGS_DEBUG} MATCHES "-Wextra")

		# fast; -fomit-frame-pointer and -finline-functions not needed, automatically set by -O3; -O4 or higher not supported by gcc;
		# -funroll-loops might be useful, but can sometimes make the code slower (see man gcc). We suppress warnings regarding unused results.
		# -O3 -Wall -std=c++98 -Wextra -pedantic -ffast-math -mtune=native -minline-stringops-dynamically -Wno-unused-result -funroll-loops -funsafe-loop-optimizations -Wunsafe-loop-optimizations"

		if(NOT ${CMAKE_CXX_FLAGS_RELEASE} MATCHES "-minline-stringops-dynamically")
			set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -std=c++98 -Wextra -pedantic -ffast-math -mtune=native -minline-stringops-dynamically -funroll-loops" CACHE "Release mode options" STRING FORCE)
		endif(NOT ${CMAKE_CXX_FLAGS_RELEASE} MATCHES "-minline-stringops-dynamically")

		#if(NOT ${CMAKE_EXE_LINKER_FLAGS_RELEASE} MATCHES "Wl")
			set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-Wl,-O1" CACHE "Release mode linker options" STRING FORCE)
		#endif(NOT ${CMAKE_EXE_LINKER_FLAGS_RELEASE} MATCHES "Wl")

	endif(CMAKE_COMPILER_IS_GNUCXX)

	if(OPENMP_FOUND)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif (OPENMP_FOUND)

	set(BLAS_LINK_DIR "/apps/eiger/Intel-FOR-11.1/mkl/lib/em64t" CACHE PATH "BLAS lib dir")
	set(LAPACK_LINK_DIR "/apps/eiger/Intel-FOR-11.1/mkl/lib/em64t" CACHE PATH "LAPACK lib dir")
	#set(MATH_LIB_NAMES blas lapack gfortranbegin gfortran CACHE STRING "Math libraries" )
	set(MATH_LIB_NAMES mkl_intel_lp64 mkl_sequential mkl_lapack mkl_core m CACHE STRING "Math libraries" )

	if(USE_LAPACK)
		add_definitions(-DUSE_LAPACK)
		link_directories( ${BLAS_LINK_DIR} )
		link_directories( ${LAPACK_LINK_DIR} )
	endif(USE_LAPACK)

	#add_definitions(-DFORCE_ARGUMENTS)

	set(NLOPT_LINK_DIR "/users/mvalle/apps/lib" CACHE PATH "NLopt lib dir")
	set(NLOPT_INCLUDE_DIR "/users/mvalle/apps/include" CACHE PATH "NLopt include dir")
	if(USE_OPTIMIZER)
		link_directories( ${NLOPT_LINK_DIR} )
		include_directories(${NLOPT_INCLUDE_DIR})
		add_definitions(-DUSE_OPTIMIZER)
	endif(USE_OPTIMIZER)

endif()

# Executable
add_executable(fast ${SRCS})

#Link libraries
if(WIN32)
	if(USE_OPTIMIZER)
		target_link_libraries(fast nlopt)
	endif(USE_OPTIMIZER)
else(WIN32)
	if(USE_OPTIMIZER)
		if(USE_LAPACK)
			target_link_libraries(fast ${MATH_LIB_NAMES} nlopt)
		else(USE_LAPACK)
			target_link_libraries(fast nlopt)
		endif(USE_LAPACK)
	else(USE_OPTIMIZER)
		if(USE_LAPACK)
			target_link_libraries(fast ${MATH_LIB_NAMES} )
		endif(USE_LAPACK)
	endif(USE_OPTIMIZER)
endif(WIN32)

# Make documentation
find_package(Doxygen)
if(DOXYGEN_FOUND)

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
   set(DOXY_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")
  endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile")

  add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG})
endif(DOXYGEN_FOUND)

