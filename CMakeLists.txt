cmake_minimum_required(VERSION 3.16)
project(parallelFrame LANGUAGES C)

# ── Options ──
option(USE_MPI "Build MPI-based targets (wave solvers, MPI examples, halo tests)" ON)

# ── Compile flags ──
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2")

# ── Dependencies ──
find_package(Threads REQUIRED)

if(USE_MPI)
    find_package(MPI REQUIRED)
    add_compile_definitions(USE_MPI)
endif()

find_package(OpenMP QUIET COMPONENTS C)

# ── mythread library ──
set(MYTHREAD_CORE
    mythread/mythread_util.c
    mythread/mythread_locv.c
    mythread/mythread_timer.c
    mythread/mythread_sync.c
    mythread/mythread_pool.c
    mythread/mythread_thread.c
    mythread/mythread_fortran.c
    mythread/mythread_decomp.c
    mythread/mythread_field.c
    mythread/mythread_config.c
)

add_library(mythread STATIC ${MYTHREAD_CORE})
target_include_directories(mythread PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(mythread PUBLIC Threads::Threads)
target_link_libraries(mythread PRIVATE m)

if(USE_MPI)
    target_sources(mythread PRIVATE mythread/mythread_halo.c)
    target_link_libraries(mythread PUBLIC MPI::MPI_C)
endif()

# ── Helper ──
macro(add_mythread_exec name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE mythread)
endmacro()

# ──────────────────────────────────────────────
# 非 MPI 目标（USE_MPI=OFF 时仍可构建）
# ──────────────────────────────────────────────
add_mythread_exec(taskpool_tests  tests/taskpool_tests.c)
add_mythread_exec(test_bindcpu    tests/test_bindcpu.c)
add_mythread_exec(test_config     tests/test_config.c)
add_mythread_exec(test_decomp     tests/test_decomp.c)
add_mythread_exec(test_locv       tests/test_locv.c)
add_mythread_exec(test_sync       tests/test_sync.c)
add_mythread_exec(test_timer      tests/test_timer.c)
add_mythread_exec(test_util       tests/test_util.c)
add_mythread_exec(example_taskpool examples/example_taskpool.c)

# ──────────────────────────────────────────────
# MPI 目标（仅 USE_MPI=ON 时构建）
# ──────────────────────────────────────────────
if(USE_MPI)
    add_mythread_exec(test_field  tests/test_field.c)
    add_mythread_exec(example_simple       examples/example_simple.c)
    add_mythread_exec(example_mpi_mythread  examples/example_mpi_mythread.c)

    add_mythread_exec(wave_propagation_ghost            wave_propagation_ghost.c)
    add_mythread_exec(wave_propagation_ghost_fix        wave_propagation_ghost_fix.c)
    add_mythread_exec(wave_propagation_ghost_mpi         wave_propagation_ghost_mpi.c)
    add_mythread_exec(wave_propagation_ghost_mpi_mythread wave_propagation_ghost_mpi_mythread.c)
    add_mythread_exec(wave_propagation_ghost_mpi_sync    wave_propagation_ghost_mpi_sync.c)

    if(OpenMP_C_FOUND)
        add_mythread_exec(wave_propagation_ghost_omp      wave_propagation_ghost_omp.c)
        target_link_libraries(wave_propagation_ghost_omp PRIVATE OpenMP::OpenMP_C)

        add_mythread_exec(wave_propagation_ghost_mpi_omp  wave_propagation_ghost_mpi_omp.c)
        target_link_libraries(wave_propagation_ghost_mpi_omp PRIVATE OpenMP::OpenMP_C)
    else()
        message(STATUS "OpenMP not found — skipping wave_propagation_ghost_omp, wave_propagation_ghost_mpi_omp")
    endif()
endif()

# ── Status ──
message(STATUS "")
message(STATUS "parallelFrame build configuration:")
message(STATUS "  MPI support:      ${USE_MPI}")
if(USE_MPI)
    message(STATUS "  MPI_CXX:          ${MPI_CXX_FOUND}")
    message(STATUS "  MPI compiler:     ${MPI_C_COMPILER}")
endif()
message(STATUS "  OpenMP:           ${OpenMP_C_FOUND}")
message(STATUS "  Threads:          ${CMAKE_THREAD_LIBS_INIT}")
message(STATUS "")