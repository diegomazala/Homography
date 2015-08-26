set (HOMOGRAPHY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Homography)

include_directories (${HOMOGRAPHY_SOURCE_DIR})

add_subdirectory (Homography)
add_subdirectory (Mosaic)
add_subdirectory (SurfOpenCV)
add_subdirectory (TestNonLinear)


