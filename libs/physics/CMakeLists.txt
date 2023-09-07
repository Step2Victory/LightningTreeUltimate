add_library(Physics physics.cpp)

target_include_directories(Physics INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Physics PUBLIC compiler_flags)

if (ENABLE_LIBS_TESTING)
    add_subdirectory(tests)
endif()