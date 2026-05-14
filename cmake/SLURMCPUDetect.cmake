# SLURM CPU Architecture Detection Module
#
# This module provides a way to automatically detect the "least common denominator"
# CPU features across all nodes in a SLURM cluster.
#
# It sets the following variables based on the common denominator:
#   NOAVX
#   NOAVX2
#   NOSSE4

option(USE_SLURM_DETECTION "Use SLURM to detect least common denominator CPU features across cluster nodes" OFF)

if(USE_SLURM_DETECTION)
    message(STATUS "SLURM CPU detection enabled. Querying cluster nodes...")
    
    find_package(Python3 REQUIRED)
    
    execute_process(
        COMMAND ${Python3_EXECUTABLE} "${CMAKE_SOURCE_DIR}/contrib/slurm_cpu_detect.py" --cmake
        OUTPUT_VARIABLE SLURM_CPU_FLAGS
        RESULT_VARIABLE SLURM_DETECTION_RESULT
        ERROR_VARIABLE SLURM_DETECTION_ERROR
    )
    
    if(SLURM_DETECTION_RESULT EQUAL 0)
        # The script outputs something like:
        # -DNOAVX2=OFF
        # -DNOAVX=OFF
        # -DNOSSE4=OFF
        
        # We need to parse these and set the corresponding CMake variables
        string(REPLACE "\n" ";" FLAG_LIST "${SLURM_CPU_FLAGS}")
        foreach(FLAG ${FLAG_LIST})
            if(FLAG MATCHES "^-D([^=]+)=(.+)$")
                set(VAR_NAME ${CMAKE_MATCH_1})
                set(VAR_VALUE ${CMAKE_MATCH_2})
                message(STATUS "SLURM Detection: Setting ${VAR_NAME}=${VAR_VALUE}")
                set(${VAR_NAME} ${VAR_VALUE} CACHE BOOL "Automatically set by SLURM detection" FORCE)
            endif()
        endforeach()
    else()
        message(WARNING "SLURM CPU detection failed: ${SLURM_DETECTION_ERROR}")
        message(STATUS "Falling back to default host-based detection.")
    endif()
endif()
