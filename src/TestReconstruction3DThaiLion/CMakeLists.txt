#######################################################################
# Setting Target_Name as current folder name
get_filename_component(TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)



set  (SOURCE_FILES		TestReconstruction3DThaiLion.cpp )

set  (HEADER_FILES						)


include_directories (${CMAKE_SOURCE_DIR})


add_executable(
  ${TARGET_NAME}
  ${SOURCE_FILES}
  ${HEADER_FILES}
)


target_link_libraries (	
	${TARGET_NAME} 
	Qt5::Widgets
	Homography
)



