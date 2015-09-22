set (HOMOGRAPHY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Homography)

include_directories (${HOMOGRAPHY_SOURCE_DIR})

add_subdirectory (Homography)
add_subdirectory (Mosaic)
add_subdirectory (Reconstruction3D)
add_subdirectory (SurfOpenCV)
add_subdirectory (TestNonLinear)

add_subdirectory (TestReconstruction3DCube)
add_subdirectory (TestReconstruction3DThaiLion)
add_subdirectory (TestReconstruction3DAlison)
add_subdirectory (TestTriangulation)
