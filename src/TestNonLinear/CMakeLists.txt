#######################################################################
# Setting Target_Name as current folder name
get_filename_component(TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)



set  (SOURCE_FILES		
						../Homography/QImageWidget.cpp 
						../Homography/DLT.cpp						
						../Homography/Points.cpp
						../Homography/RansacDLT.cpp
						../Homography/GaussNewton.cpp
						TestNonLinear.cpp )

set  (HEADER_FILES		
						../Homography/QImageWidget.h
						../Homography/ImageHelper.h
						../Homography/DLT.h
						../Homography/Points.h
						../Homography/RansacDLT.h
						../Homography/GaussNewton.h
						)

#set  (FORM_FILES		 
#						MainWindow.ui)


#qt5_wrap_ui( HEADER_UI_FILES ${FORM_FILES} )

source_group ("Form Files" FILES ${FORM_FILES})
source_group ("Form Header Files" FILES ${HEADER_UI_FILES})


include_directories (${CMAKE_SOURCE_DIR} ../Homography/)


add_executable(
  ${TARGET_NAME}
  ${SOURCE_FILES}
  ${HEADER_FILES}
  ${HEADER_UI_FILES}
)

target_link_libraries (	
	${TARGET_NAME} 
	Qt5::Widgets
)



