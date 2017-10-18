find_path(LZ4_INCLUDE_DIR lz4.h HINTS ${LZ4_ROOT})

find_library(LZ4_LIBRARY NAMES lz4 HINTS ${LZ4_ROOT})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LZ4 DEFAULT_MSG LZ4_LIBRARY LZ4_INCLUDE_DIR)

mark_as_advanced(LZ4_INCLUDE_DIR LZ4_LIBRARY)

set(LZ4_LIBRARIES ${LZ4_LIBRARY})
set(LZ4_INCLUDE_DIRS ${LZ4_INCLUDE_DIR})
