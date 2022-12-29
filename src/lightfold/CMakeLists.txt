cmake_minimum_required(VERSION 3.13)

add_executable(lightfold lightfold.cpp)

target_include_directories(
    lightfold
    PUBLIC ${TINYEXR_INCLUDE_DIR}
)