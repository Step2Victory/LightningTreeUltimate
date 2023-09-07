add_library(Libs INTERFACE)

add_subdirectory(physics)

target_link_libraries(Libs INTERFACE Physics)