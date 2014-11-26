#####
# WINDOWS USERS: The easiest way to make cmake find boost/Eigen/OpenCV is to provide 
# all of the paths to library and include dirs here.
# 
# If anyone knows how to make this work please let me know :)
#
#MESSAGE(STATUS "CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}")
IF(WIN32) # (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  MESSAGE(STATUS "The easiest way to make cmake find boost/Eigen/OpenCV is to provide all of the paths to library and include dirs here (NB use FORWARD SLASHES).")

  set(Boost_INCLUDE_DIR "G:/Users/Tom Botterill/Documents/lib/boost_1_55_0")
  set(OpenCV_INCLUDE_DIR "G:/Users/Tom Botterill/Documents/lib/opencv/opencv/build/include")
  set(EIGEN3_INCLUDE_DIR "G:/Users/Tom Botterill/Documents/lib/eigen")
  set(Boost_LIBRARY_DIRS  "G:/Users/Tom Botterill/Documents/lib/boost_1_55_0/stage/lib")
  set(Boost_FOUND  TRUE)
  set(OpenCV_LIBRARIES "G:/Users/Tom Botterill/Documents/lib/opencv/opencv/build/x64/vc12/lib")
endif()

######################


cmake_minimum_required (VERSION 2.6)
project(tom-cv)

include_directories(util)
include_directories(timing)
include_directories(params)
include_directories(logging)
include_directories(cameraGeom)
include_directories(ransac)
include_directories(bagOfWords)
include_directories(bowSpeedometer)
include_directories(featureDescription)
include_directories(featureExtract)
include_directories(featureMatching)
include_directories(image)
include_directories(imageSource)


#################################################### check for boost, Eigen3, opencv

 # Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")

  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
  endif()    
endif()

set(LIBRARY_OUTPUT_PATH "../${CMAKE_BUILD_TYPE}")


IF(NOT Boost_FOUND)
  find_package (Boost 1.54.0 REQUIRED COMPONENTS system filesystem thread)
ENDIF()

include_directories (SYSTEM ${Boost_INCLUDE_DIR})
message(STATUS "boost include dir= ${Boost_INCLUDE_DIR}")

#############

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ./cmake)

IF(NOT OpenCV_INCLUDE_DIR)
  find_package(OpenCV 2.4 REQUIRED COMPONENTS opencv_core opencv_highgui opencv_imgproc opencv_legacy opencv_ml opencv_calib3d) 
endif()

message(STATUS "opencv include dir= ${OpenCV_INCLUDE_DIR}")
message(STATUS "opencv root dir+include= ${OpenCV_ROOT_DIR}/include")


message(STATUS "opencv_found= ${OpenCV_FOUND}")

include_directories(SYSTEM ${OpenCV_INCLUDE_DIR})
include_directories(SYSTEM "${OpenCV_ROOT_DIR}/include")

#############
IF(NOT EIGEN3_INCLUDE_DIR)
  FIND_PACKAGE(Eigen3 3.1 REQUIRED)
ENDIF()

message(STATUS "Eigen3 include dir= ${EIGEN3_INCLUDE_DIR}")

message(STATUS "Eigen3_found= ${EIGEN3_FOUND}")

include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

#################################################### build all the tom-cv libraries

add_subdirectory(util)
add_subdirectory(timing)
add_subdirectory(params)
add_subdirectory(logging)
add_subdirectory(cameraGeom)
add_subdirectory(ransac)
add_subdirectory(bagOfWords)
add_subdirectory(bowSpeedometer)
add_subdirectory(featureDescription)
add_subdirectory(featureExtract)
add_subdirectory(featureMatching)
add_subdirectory(image)
add_subdirectory(imageSource)

################################################### set up libraries to link to
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")

  link_libraries(featureExtract)

  link_libraries(imageSource)

  link_libraries(params)
  link_libraries(logging)
  link_libraries(ransac)
  link_libraries(bagOfWords)
  link_libraries(bowSpeedometer)
  link_libraries(featureDescription)
  link_libraries(featureMatching)
  link_libraries(image)

  link_libraries(${Boost_LIBRARIES})
  link_libraries(opencv_core)
  link_libraries(opencv_imgproc)
  link_libraries(opencv_calib3d)
  link_libraries(opencv_highgui)
  link_libraries(opencv_features2d)
  link_libraries(opencv_legacy)


  link_libraries(cameraGeom)
  link_libraries(timing)



  link_libraries(util)

  link_libraries(png)
  link_libraries(jpeg)


################################################### build the example programs

  add_subdirectory(ransacTest)
  add_subdirectory(ransacAndRefine)
  add_subdirectory(bowExample)

  add_subdirectory(mosaicing)
  add_subdirectory(BoWSLAM)
  add_subdirectory(ShowIm) # needed for BoWSLAM

ENDIF()
