# MIT License
#
# Copyright (c) 2019 University of Oxford
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.6)

project(
    DistributionZoo
    VERSION 0.1
    DESCRIPTION "A zoo of univariate distributions."
    LANGUAGES CXX
)

# ensure C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#export compiler flags for code completion engines
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Set global compiler warnings
if (MSVC)
    add_compile_options(/W1 /WX)
else ()
    add_compile_options(-Wall -pedantic)
endif ()

add_library(cts_univ INTERFACE)
target_include_directories(cts_univ INTERFACE continuous_univariate)

add_library(dsc_univ INTERFACE)
target_include_directories(dsc_univ INTERFACE discrete_univariate)

add_library(zoo_util INTERFACE)
target_include_directories(zoo_util INTERFACE zoo_util)

enable_testing()


set(
        TEST_FILES
        tests/tests_main.cpp
        tests/continuous_univariate_tests.cpp
        tests/discrete_univariate_tests.cpp
        tests/zoo_util_tests.cpp
)

add_executable(tests ${TEST_FILES})
target_link_libraries(tests PRIVATE cts_univ)
target_link_libraries(tests PRIVATE dsc_univ)
target_link_libraries(tests PRIVATE zoo_util)
add_test(tests tests)

# Clang tidy as optional static analyzer
option(Zoo_USE_CLANG_TIDY "Use clang tidy for static analysis" OFF)
if (Zoo_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-7" "clang-tidy-6.0"
            DOC "Path to clang-tidy executable")
    if (NOT CLANG_TIDY_EXE)
        message(STATUS "clang-tidy not found.")
    else ()
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set_target_properties(cts_univ PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    endif ()
endif ()

# LLVM AddressSanitizer to detecting memory errors
# Note that there are many other sanitizers in LLVM to help detect errors, see
# http://travistoptips.blogspot.co.uk/2015/11/sanitize-all-things.html
option(Zoo_MEMCHECK "Use LLVM AddressSanitizer to detecting memory errors" OFF)
if (Zoo_MEMCHECK)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        message(STATUS "Configuring with LLVM AddressSanitizer")
        set(Zoo_MEMCHECK_FLAGS -fno-optimize-sibling-calls
                -fsanitize=address
                -fsanitize-address-use-after-scope
                )
        target_compile_options(cts_univ INTERFACE -O1 -g -fno-omit-frame-pointer ${Zoo_MEMCHECK_FLAGS})
        target_link_libraries(cts_univ INTERFACE -g ${Zoo_MEMCHECK_FLAGS})

        target_compile_options(dsc_univ INTERFACE -O1 -g -fno-omit-frame-pointer ${Zoo_MEMCHECK_FLAGS})
        target_link_libraries(dsc_univ INTERFACE -g ${Zoo_MEMCHECK_FLAGS})

        target_compile_options(zoo_util INTERFACE -O1 -g -fno-omit-frame-pointer ${Zoo_MEMCHECK_FLAGS})
        target_link_libraries(zoo_util INTERFACE -g ${Zoo_MEMCHECK_FLAGS})
    else ()
        message(FATAL_ERROR "clang compiler required with Zoo_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()

# Setup coverage testing for GCC or Clang
option(Zoo_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (Zoo_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message(STATUS "Configuring with coverage")
        target_compile_options(cts_univ INTERFACE --coverage -O0)
        target_link_libraries(cts_univ INTERFACE --coverage)

        target_compile_options(dsc_univ INTERFACE --coverage -O0)
        target_link_libraries(dsc_univ INTERFACE --coverage)

        target_compile_options(zoo_util INTERFACE --coverage -O0)
        target_link_libraries(zoo_util INTERFACE --coverage)
    else ()
        message(FATAL_ERROR "GCC or Clang required with Zoo_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()
