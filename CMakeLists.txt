cmake_minimum_required(VERSION 3.13)
project(graph_balls_cpp)

set(CMAKE_CXX_STANDARD 17)

foreach (OUTPUTCONFIG IN ITEMS Debug Release)
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_SOURCE_DIR}/build-${OUTPUTCONFIG}/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_SOURCE_DIR}/build-${OUTPUTCONFIG}/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_SOURCE_DIR}/build-${OUTPUTCONFIG}/lib")
endforeach ()

option(USE_TBBMALLOC "Use TBB memory allocator" OFF)

set(EXECUTABLE_NAME graphballs)

# bliss creator decided that they're too cool for && and ||
if (MSVC)
    # /permissive- for standard C++ with "and" and "or" logical operators
    add_compile_options(/permissive-)
else()
    # Warnings, optimization, no assertions
    add_compile_options(-Wall -pedantic -O3 -DNDEBUG) #-Wextra -Werror
endif()

include_directories(src include)
link_directories(lib)

add_executable(${EXECUTABLE_NAME}
        src/main.cpp
        src/GraphAdjacency.h
        src/GraphAdjacency.cpp
        src/GraphBallIsomorphism.cpp
        src/GraphBallIsomorphism.h
        src/GraphPartitioning.cpp
        src/GraphPartitioning.h
        src/GraphBallsUtil.cpp
        src/GraphBallsUtil.h
        src/GraphPartitioningThreadPool.cpp
        src/GraphPartitioningThreadPool.h
)
add_subdirectory(external/bliss-0.77)
target_include_directories(${EXECUTABLE_NAME} PRIVATE external/bliss-0.77/src)
target_link_libraries(${EXECUTABLE_NAME} PRIVATE bliss_static)

# TBB malloc replaces memory system calls with a custom allocator which makes allocations faster, especially multithread
if (USE_TBBMALLOC)
    message("use tbbmalloc")
    if (MSVC)
        SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} tbbmalloc_proxy.lib /INCLUDE:\"__TBB_malloc_proxy\"")
    endif ()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        add_library(tbbmalloc_proxy SHARED IMPORTED)
        set_target_properties(tbbmalloc_proxy PROPERTIES IMPORTED_LOCATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/libtbbmalloc_proxy.so)
        target_link_libraries(${EXECUTABLE_NAME} tbbmalloc_proxy)
    endif ()
else ()
    message("not use tbbmalloc")
endif ()
unset(USE_TBBMALLOC CACHE)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE pthread stdc++fs)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif ()
