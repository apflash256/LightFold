cmake_minimum_required(VERSION 3.12)
project(LightFold LANGUAGES C CXX)

# ext
add_subdirectory(src/ext)

# LensSim - include and link

#target_include_directories(lightfold
#  PUBLIC
#    src
#    src/lightfold
#)
#target_link_libraries(lightfold PUBLIC tinyexr)

add_executable(lightfold src/lightfold/lightfold.cpp)
set_target_properties(lightfold  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(lightfold PRIVATE ${CMAKE_SOURCE_DIR}/src/ext)
target_include_directories(lightfold PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(lightfold PRIVATE tinyexr::tinyexr)