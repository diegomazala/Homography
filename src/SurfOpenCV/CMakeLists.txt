
project( surf_matcher )
add_executable( surf_matcher surf_matcher.cpp )
target_link_libraries( surf_matcher ${OpenCV_LIBS} )