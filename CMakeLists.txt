cmake_minimum_required(VERSION 3.0...3.8)

####
# project settings
####

project(DemoObs)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	main.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(DemoObs ${HEADERS} ${SOURCES})

target_link_libraries(DemoObs ${GLAMER_LIBRARIES})


