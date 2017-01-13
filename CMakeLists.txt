# Some functions require a certain version of cmake (for example file GLOB )
cmake_minimum_required(VERSION 2.8)

# Build options
option(BUILD_OPENGL_EXAMPLE "Build OpenGL example." OFF)

# Create the project
set(PROJECT_NAME fluidsim)
project(${PROJECT_NAME})

# Add the external module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

# Find external packages
find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories(${EIGEN3_INCLUDE_DIRS})

# Our own include directory
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${COMMON_INCLUDES})

# Debugger flag
set(CMAKE_CXX_FLAGS "-g")

# Creating our own ${PROJECCT_NAME} library
file(GLOB ${PROJECCT_NAME}_lib_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB ${PROJECCT_NAME}_lib_HEADERS ${PROJECT_SOURCE_DIR}/include/*.h)
add_library(${PROJECT_NAME}_lib
	${${PROJECCT_NAME}_lib_SOURCES}
	${${PROJECCT_NAME}_lib_HEADERS})
set_target_properties(${PROJECT_NAME}_lib PROPERTIES COMPILE_FLAGS "-std=c++14")

# Add program executable (all examples)
# Simple example
add_executable(${PROJECT_NAME}_simple examples/simple.cpp)
target_link_libraries(${PROJECT_NAME}_simple
	${PROJECT_NAME}_lib
	)
set_target_properties(${PROJECT_NAME}_simple PROPERTIES COMPILE_FLAGS "-std=c++14")

# OpenGL example
if(BUILD_OPENGL_EXAMPLE)
	find_package(GLFW3 REQUIRED)
  	include_directories(${GLFW_INCLUDE_DIRS})

	add_subdirectory(${PROJECT_SOURCE_DIR}/ext/Simple_Graphics_Engine)
	include_directories(${SGE_INCLUDE_DIRS})

	file(GLOB ${PROJECCT_NAME}_openGL_example_SOURCES
		${PROJECT_SOURCE_DIR}/examples/opengl_example/*.cpp)
	file(GLOB ${PROJECCT_NAME}_openGL_example_HEADERS
		${PROJECT_SOURCE_DIR}/examples/opengl_example/*.h)

	add_executable(${PROJECT_NAME}_openGL_example
		${${PROJECCT_NAME}_openGL_example_SOURCES}
		${${PROJECCT_NAME}_openGL_example_HEADERS})
	target_link_libraries(${PROJECT_NAME}_openGL_example
		${PROJECT_NAME}_lib
		SGE
		${GLFW_LIBRARIES}
		)
	set_target_properties(${PROJECT_NAME}_openGL_example
		PROPERTIES COMPILE_FLAGS "-std=c++14")
	target_compile_definitions(${PROJECT_NAME}_openGL_example
		PRIVATE PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")

endif()