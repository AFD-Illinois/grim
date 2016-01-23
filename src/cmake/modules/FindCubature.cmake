# - Try to find Cubature source code and headers.
#
# Usage of this module as follows:
#
#     find_package(Cubature)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  CUBATURE_DIR  Set this variable to the root installation of
#                CUBATURE if the module has problems finding
#                the proper installation path.
#
# Variables defined by this module:
#
# CUBATURE_FOUND           System has Cubature headers
# CUBATURE_INCLUDE_DIR     The location of Cubature headers

find_path(CUBATURE_DIR
    NAMES cubature.h
)

find_path(CUBATURE_INCLUDE_DIR
    NAMES cubature.h
    HINTS ${CUBATURE_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CUBATURE DEFAULT_MSG
    CUBATURE_INCLUDE_DIR
)

mark_as_advanced(
    CUBATURE_DIR
    CUBATURE_INCLUDE_DIR
)
