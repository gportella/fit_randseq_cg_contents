# Minimum cmake version
cmake_minimum_required (VERSION 3.5.0)

# Name of project and that it is C++ only.
set(code_name "cg_content")
project (${code_name} )

# ----------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------- 
# Search for zlib as a dependency for SeqAn.
set(CMAKE_MODULE_PATH ~/work/soft/seqan-library-2.2.0/share/cmake/Modules) 
set(CMAKE_PREFIX_PATH ~/work/soft/seqan-library-2.2.0)
find_package (ZLIB)
find_package (BZip2)
# Would be nice to add an ifdef depending on the machie
execute_process(COMMAND hostname OUTPUT_VARIABLE computer_name)
find_package (SeqAn REQUIRED)



# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS}  )

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).

# Do use -O3, makes the code way faster
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g  ")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -g -Wall -Wextra -O3 ")

add_executable (${code_name} ${code_name}.cpp utils_common.hpp compute_contents.hpp  )
target_link_libraries (${code_name} ${SEQAN_LIBRARIES})

