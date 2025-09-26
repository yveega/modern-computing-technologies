cmake_minimum_required(VERSION 3.5)
project(inmost_project)

#set(INMOST_DIR %%%%%) # uncomment and set location manually if needed
find_package(inmost REQUIRED)
find_package(Matplot++ REQUIRED)

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(main main.cpp)
target_include_directories(main PUBLIC src)
target_link_libraries(main PUBLIC ${INMOST_LIBRARIES})

add_executable(visualize visualize.cpp)
target_include_directories(visualize PUBLIC src)
target_link_libraries(visualize PUBLIC Matplot++::matplot)
target_link_libraries(visualize PUBLIC ${INMOST_LIBRARIES})