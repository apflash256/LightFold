cmake_minimum_required(VERSION 3.12)
project(LightFold LANGUAGES C CXX)

# ext
add_subdirectory(src/ext)

# src
add_subdirectory(src/lightfold)