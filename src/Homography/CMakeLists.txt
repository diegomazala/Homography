#######################################################################
# Setting Target_Name as current folder name
get_filename_component(TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)



set  (SOURCE_FILES		
						QImageWidget.cpp 
						DLT.cpp						
						Points.cpp
						RansacDLT.cpp
						AppMain.cpp )

set  (HEADER_FILES		
						QImageWidget.h
						ImageHelper.h
						DLT.h
						Points.h
						RansacDLT.h
						)

#set  (FORM_FILES		 
#						MainWindow.ui)


#qt5_wrap_ui( HEADER_UI_FILES ${FORM_FILES} )

source_group ("Form Files" FILES ${FORM_FILES})
source_group ("Form Header Files" FILES ${HEADER_UI_FILES})


include_directories (${CMAKE_SOURCE_DIR})


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



