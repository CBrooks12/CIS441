cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(finalProject)
SET(VTK_DIR /path/to/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(finalProject finalProject.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(finalProject ${VTK_LIBRARIES})
else()
  target_link_libraries(finalProject vtkHybrid)
endif()


