cmake_minimum_required(VERSION 3.13)

project(lightfold)

### Debug
#message("CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
#message("CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
#message("CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})
#message("CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR})

### External
add_subdirectory(src/ext)
add_subdirectory(src/lightfold)