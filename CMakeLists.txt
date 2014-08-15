project(WaveguideModeSolver)
	cmake_minimum_required(VERSION 2.8)

	set(CMAKE_MODULE_PATH ${CMAKE_HOME_DIRECTORY}/CMakeMacros)
	
	# Allow to override the default values in CMake.vars:
	include(CMake.vars OPTIONAL)
      
    find_package(HERMES REQUIRED)
	include_directories(${HERMES_INCLUDE})

  	find_package(HERMES_COMMON REQUIRED)
	include_directories(${HERMES_COMMON_INCLUDE})

	find_package(Eigen3 3.2 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})
	add_definitions(-DEIGEN_USE_MKL_ALL)

	#find_package(MKL REQUIRED)

	find_package(Boost 1.36.1 REQUIRED COMPONENTS program_options)
	IF (Boost_FOUND)
		INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
		ADD_DEFINITIONS( "-DHAS_BOOST" )
		MESSAGE (STATUS "BOOST includes found in " ${Boost_INCLUDE_DIR})
		MESSAGE (STATUS "BOOST libraries found in " ${Boost_LIBRARY_DIRS})
	ENDIF()

	#find_package(ARPACK REQUIRED)

	#find_package(UMFPACK REQUIRED)
	#include_directories(${UMFPACK_INCLUDE_DIR})

	## === Now, some configuration for the various compilers ===

	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
			# Disable all warnings and turn on only important ones:
			set(CMAKE_CXX_FLAGS " -std=c++0x -Wno-unknown-pragmas ${CMAKE_CXX_FLAGS}")
			set(RELEASE_FLAGS "-DNDEBUG -O3")
			set(DEBUG_FLAGS   "-g")
	else()
			# Disable all warnings and turn on only important ones:
			set(CMAKE_CXX_FLAGS " -std=c++11 -Wno-unknown-pragmas ${CMAKE_CXX_FLAGS}")
			set(RELEASE_FLAGS "-DNDEBUG -Ofast")
			set(DEBUG_FLAGS   "-g")
	endif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

	# Enabling multiprocessor build on MSVC
	if(MSVC)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") 
	endif(MSVC)

	add_subdirectory(Source)

  
  
  
