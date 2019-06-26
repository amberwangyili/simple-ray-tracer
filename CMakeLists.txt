cmake_minimum_required(VERSION 2.8)
project(tracer_final)

find_package(OpenCV)
find_package(OpenMP)
add_definitions(-std=c++11)
add_definitions(-O2)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories( ${OpenMP_INCLUDE_DIRS} )


set(CMAKE_CXX_STANDARD 14)

add_executable(tracer_final main.cpp Config.h Scene.h tracer/Pathtracer.h tracer/PPMtracer.h surface/Bezier.h surface/Implicit.h kdtree/Kdtree.h object/Material.h object/Object.h object/Texture.h base/Camera.h base/Canvas.h base/Common.h base/Random.h base/Ray.h base/Vec.h kdtree/PhotonKdtree.h tracer/Photon.h)

target_link_libraries(  tracer_final ${OpenCV_LIBS} ${OpenMP_LIBS})
