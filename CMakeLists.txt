cmake_minimum_required(VERSION 3.5)

project(integral ${PROJECT_NAME} LANGUAGES CXX )

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Boost ${Boost_VERSION_MAJOR} COMPONENTS unit_test_framework REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file(GLOB SOURCES
    solver.hpp
    main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})


target_link_libraries (${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX  ${Boost_LIBRARIES})

file(GLOB SOURCES_
    solver.hpp
    )

set (TESTS_SOURCES unit_test.cpp)
set (TEST test_${PROJECT_NAME})
add_executable (${TEST} ${TESTS_SOURCES} ${SOURCES_})
target_link_libraries (${TEST}  PUBLIC OpenMP::OpenMP_CXX ${Boost_LIBRARIES})
enable_testing ()
add_test (${TEST} ${TEST})
