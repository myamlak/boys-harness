cmake_minimum_required(VERSION 3.20)

project(d_second_instance LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# The claimgate package this instance is built against is installed from the
# v0.1.0 tag; find_package reads that package's own version file.
find_package(claimgate 0.1.0 CONFIG REQUIRED)

add_executable(instance src/instance.cpp)
target_link_libraries(instance PRIVATE claimgate::claimgate)

if(MSVC)
    target_compile_options(instance PRIVATE /W4 /permissive-)
else()
    target_compile_options(instance PRIVATE -Wall -Wextra -Wpedantic)
endif()
