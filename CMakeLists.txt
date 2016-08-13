cmake_minimum_required(VERSION 3.0)
project(2d-lbm-dem C)
set(EXECUTABLE_OUTPUT_PATH "bin")
add_executable(lbmdem src/main.c)
target_link_libraries(lbmdem m)

option (USE_AUTOPARALLEL  "Use compiler auto-parallelisation" ON)
option (USE_OPENMP        "Use OpenMP directives" ON)
option (USE_OPENACC       "Use OpenACC directives" ON)


set (C_FLAGS_GNU   "-std=c99 -g")
set (C_FLAGS_Intel "-std=c99 -g")
set (C_FLAGS_PGI   "-std=c99 -g")
set (C_FLAGS "${C_FLAGS_${CMAKE_C_COMPILER_ID}}")

set (C_FLAGS_RELEASE_GNU   "-Ofast -lm -std=c99 -march=native -flto -fipa-pta -funsafe-math-optimizations -fno-math-errno -fno-trapping-math")
set (C_FLAGS_RELEASE_Intel "-O3 -xHost -ipo -fp-model fast=2 -no-prec-div")
set (C_FLAGS_RELEASE_PGI   "-O3 -Mipa")
set (C_FLAGS_RELEASE "${C_FLAGS_RELEASE_${CMAKE_C_COMPILER_ID}}")

set (C_FLAGS_DEBUG_GNU   "-Og -Wall -Wextra")
set (C_FLAGS_DEBUG_Intel "")
set (C_FLAGS_DEBUG_PGI   "")
set (C_FLAGS_DEBUG "${C_FLAGS_DEBUG_${CMAKE_C_COMPILER_ID}}")

set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${C_FLAGS} ${C_FLAGS_RELEASE}")
set (CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG}   ${C_FLAGS} ${C_FLAGS_DEBUG}")


set (C_FLAGS_AUTOPARALLEL_Intel "-parallel")
set (C_FLAGS_AUTOPARALLEL_PGI   "-Mconcur")
set (C_FLAGS_AUTOPARALLEL "${C_FLAGS_AUTOPARALLEL_${CMAKE_C_COMPILER_ID}}")

set (C_FLAGS_OPENMP_GNU   "-fopenmp")
set (C_FLAGS_OPENMP_Intel "-qopenmp")
set (C_FLAGS_OPENMP_PGI   "-mp")
set (C_FLAGS_OPENMP "${C_FLAGS_OPENMP_${CMAKE_C_COMPILER_ID}}")

set (C_FLAGS_OPENACC_PGI "-acc")
set (C_FLAGS_OPENACC "${C_FLAGS_OPENACC_${CMAKE_C_COMPILER_ID}}")

if (USE_AUTOPARALLEL)
  message("== Building with auto-parallelisation")
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${C_FLAGS_AUTOPARALLEL}")
endif()
if (USE_OPENMP)
  message("== Building with OpenMP")
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${C_FLAGS_OPENMP}")
endif()
if (USE_OPENACC)
  message("== Building with OpenACC")
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${C_FLAGS_OPENACC}")
endif()

