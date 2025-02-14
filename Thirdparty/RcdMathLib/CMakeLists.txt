# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(RcdMathLib VERSION 1.0 LANGUAGES C)

# Set the C standard to C11
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_C_EXTENSIONS On)

set(ProjectType STATIC) # SHARED or STATIC

# Add the library
set(SourceFiles
    "src/combinatorics.c"
    "src/damped_newton_raphson.c"
    "src/fsolve.c"
    "src/levenberg_marquardt.c"
    "src/lu_decomp.c"
    "src/matrix.c"
    "src/modified_gauss_newton.c"
    "src/moore_penrose_pseudo_inverse.c"
    "src/newton_raphson.c"
    "src/norm_dist_rnd_generator.c"
    "src/qr_common.c"
    "src/qr_givens.c"
    "src/qr_householder.c"
    "src/qr_pseudo_inverse.c"
    "src/shell_sort.c"
    "src/solve.c"
    "src/svd.c"
    "src/utils.c"
    "src/vector.c"
)

include_directories(include)
add_library(${PROJECT_NAME} ${ProjectType} ${SourceFiles})

# Create alias target for consistent naming
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(RcdMathLib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(GNUInstallDirs)
install(TARGETS RcdMathLib
    EXPORT RcdMathLibTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install header files with rcd/ prefix
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rcd
    FILES_MATCHING PATTERN "*.h"
)

install(EXPORT RcdMathLibTargets
    FILE RcdMathLibTargets.cmake
    NAMESPACE RcdMathLib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RcdMathLib
)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/RcdMathLibConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RcdMathLib
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/RcdMathLibConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RcdMathLib
)