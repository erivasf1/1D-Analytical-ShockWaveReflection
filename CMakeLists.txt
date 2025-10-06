cmake_minimum_required(VERSION 3.10)

# Project name
project(1DFSRP)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_PREFIX_PATH "/apps/easybuild/software/owl-genoa/Catch2/2.13.10-GCCcore-13.3.0/lib")
# Fetch Catch2 in the system
find_package(Catch2 REQUIRED)


# Set the default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -lmpi")
##set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -fsanitize=address -lmpi")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



# Define source files
set(SOURCES
    Main.cpp
    Visualizer.cpp
    RankineHugoniot.cpp
)


# Define header files
set(HEADERS
    Visualizer.h
    RankineHugoniot.h
)


# Create the executable
add_executable(fsrp ${SOURCES} ${HEADERS})

# Include directories (if needed)
#target_include_directories(NozzleSimulation PRIVATE ${CMAKE_SOURCE_DIR})

# Link libraries (if any)
#target_link_libraries(nozzle PRIVATE ${CMAKE_PREFIX_PATH})
#target_link_libraries(NozzleSimulation <LIBRARY_NAME>)
