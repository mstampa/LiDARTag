cmake_minimum_required (VERSION 3.1)
project (lidartag)

# CMAKE TWEAKS
#========================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -DEIGEN_NO_DEBUG -march=native -Wl,--no-as-needed")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
message(STATUS "=============================================CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    tf
    velodyne_pcl
    velodyne_pointcloud
    visualization_msgs
    lidartag_msgs
)

# PCL
find_package(PCL 1.2 REQUIRED)
find_package(Boost REQUIRED
  COMPONENTS filesystem system regex date_time program_options thread
)

# NLOPT
find_package(NLopt REQUIRED)
message(STATUS "===============================NLOPT lib: ${NLOPT_LIBRARIES}")

# TBB
find_package(TBB REQUIRED)
if(TBB_FOUND)
    message(STATUS "=============================================TBB FOUND")
    message(STATUS "===============================TBB include_dirs: ${TBB_INCLUDE_DIRS}")
    message(STATUS "===============================TBB includes: ${TBB_INCLUDES}")
    message(STATUS "===============================TBB libs: ${TBB_LIBS}")
    message(STATUS "===============================TBB libraries: ${TBB_LIBRARIES}")
    message(STATUS "===============================TBB defs: ${TBB_DEFINITIONS}")
else()
    message(STATUS "=============================================TBB NOT FOUND")
endif() 

message(STATUS "Trying to include Eigen library")
set(EIGEN3_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/internal_eigen3/eigen)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Eigen3 3.4.0 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message(STATUS "================================== EIGEN (internal): " ${EIGEN3_VERSION})

catkin_package(
  INCLUDE_DIRS include
  CATKIN_DEPENDS 
    roscpp
    tf
    velodyne_pcl
    velodyne_pointcloud
    visualization_msgs
    lidartag_msgs
)

include_directories(
	include 
	${catkin_INCLUDE_DIRS} 
	${PCL_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${TBB_INCLUDE_DIRS}
)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# COMPILE THE SOURCE
#========================================================================
add_executable(lidartag_main
  src/apriltag_utils.cc
  src/lidartag.cc 
  src/lidartag_cluster.cc 
  src/lidartag_decode.cc
  src/lidartag_pose.cc 
  src/lidartag_prune.cc
  src/lidartag_rviz.cc
  src/main.cc 
  src/tag49h14.cc
  src/tag16h5.cc
  src/utils.cc
)

add_dependencies(lidartag_main 
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(lidartag_main 
  ${catkin_LIBRARIES} 
  ${Boost_LIBRARIES}
  ${NLOPT_LIBRARIES}
  ${PCL_LIBRARIES}
  ${TBB_LIBRARIES}
)
