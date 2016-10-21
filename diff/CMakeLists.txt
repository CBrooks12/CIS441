cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(differencer)
SET(VTK_DIR VTK)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(differencer differencer)

if(VTK_LIBRARIES)
  target_link_libraries(differencer ${VTK_LIBRARIES})
else()
  target_link_libraries(differencer vtkHybrid)
endif()

