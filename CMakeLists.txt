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
project(DistributionZoo)

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

add_library(univariate src/univariate_distributions.hpp)
target_include_directories(univariate PUBLIC src)

set_target_properties(univariate PROPERTIES LINKER_LANGUAGE CXX)


enable_testing()

add_executable(univ_tests tests/tests.cpp)
target_link_libraries(univ_tests PRIVATE univariate)
add_test(univ_test univ_tests)

# Clang tidy as optional static analyzer
option(Zoo_USE_CLANG_TIDY "Use clang tidy for static analysis" OFF)
if (Zoo_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy" "clang-tidy-7" "clang-tidy-6.0"
            DOC "Path to clang-tidy executable")
    if (NOT CLANG_TIDY_EXE)
        message(STATUS "clang-tidy not found.")
    else ()
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set_target_properties(univariate PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
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
        target_compile_options(univariate PUBLIC -O1 -g -fno-omit-frame-pointer ${Zoo_MEMCHECK_FLAGS})
        target_link_libraries(univariate PUBLIC -g ${Zoo_MEMCHECK_FLAGS})
    else ()
        message(FATAL_ERROR "clang compiler required with Zoo_MEMCHECK: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()

# Setup coverage testing for GCC or Clang
option(Zoo_ENABLE_COVERAGE "Enable coverage reporting for GCC or Clang" FALSE)
if (Zoo_ENABLE_COVERAGE)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        message(STATUS "Configuring with coverage")
        target_compile_options(univariate PUBLIC --coverage -O0)
        target_link_libraries(univariate PUBLIC --coverage)
    else ()
        message(FATAL_ERROR "GCC or Clang required with Zoo_ENABLE_COVERAGE: found ${CMAKE_CXX_COMPILER_ID}")
    endif ()
endif ()
