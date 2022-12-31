file(GLOB_RECURSE LIGHTFOLD_SRC "*.cpp")

add_executable(lightfold ${LIGHTFOLD_SRC})
set_target_properties(lightfold  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(lightfold
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/ext
  PUBLIC
    ${CMAKE_SOURCE_DIR}/src/lightfold
)
target_link_libraries(lightfold PRIVATE tinyexr::tinyexr)