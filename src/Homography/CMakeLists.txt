#######################################################################
# Setting Target_Name as current folder name
get_filename_component(TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)



set  (SOURCE_FILES		
						QImageWidget.cpp 
						DLT.cpp						
						Points.cpp
						RansacDLT.cpp
						GaussNewton.cpp
						Reconstruction3D.cpp 
						Triangulation.cpp)

set  (HEADER_FILES		
						QImageWidget.h
						ImageHelper.h
						DLT.h
						Points.h
						GaussNewton.h
						RansacDLT.h
						Reconstruction3D.h
						Triangulation.h
						ObjHelper.h)


include_directories (${CMAKE_SOURCE_DIR})


add_library(
  ${TARGET_NAME}
  ${SOURCE_FILES}
  ${HEADER_FILES}
)


target_link_libraries (	
	${TARGET_NAME} 
)



