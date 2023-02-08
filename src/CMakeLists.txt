add_executable (cutStock "main.cpp")
find_package(ortools REQUIRED CONFIG)
target_link_libraries(cutStock PRIVATE ortools::ortools)
