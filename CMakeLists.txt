cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)

project(LightningTreeUltimate)

find_package (glog 0.6.0 REQUIRED)
find_package (yaml-cpp REQUIRED)

add_executable(main src/LightningTree.cpp main.cpp)
target_link_libraries (main glog::glog)
target_link_libraries (main yaml-cpp)