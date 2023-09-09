cmake_minimum_required(VERSION 3.15)

project(LightningTreeUltimate)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-O3;-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
  "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

option(ENABLE_LIBS_TESTING "Включить модульное тестирование" ON)
message(STATUS "testing " ${ENABLE_LIBS_TESTING})
add_subdirectory(libs)

find_package(glog 0.6.0 REQUIRED)
find_package(yaml-cpp REQUIRED)


add_executable(main main.cpp src/LightningTree.cpp)
target_link_libraries(main glog::glog
                            yaml-cpp
                            compiler_flags)

add_executable(LibsExample libs_example.cpp)
target_link_libraries(LibsExample Libs)