PROJECT(ImageProcessing)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(CMAKE_CXX_FLAGS "-fopenmp -O0 -Wall -fPIC -std=c++0x ${CMAKE_CXX_FLAGS}")
#SET(CMAKE_BUILD_TYPE "Release")

ADD_SUBDIRECTORY(CommonFile)

# EIGEN3
FIND_PACKAGE(Eigen3 QUIET REQUIRED)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
  MESSAGE(STATUS "Found EIGEN3 @ ${EIGEN3_INCLUDE_DIR}")
ELSE(EIGEN3_FOUND)
  MESSAGE(SEND_ERROR "Cannot find eigen3")
ENDIF(EIGEN3_FOUND)

# BOOST
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost QUIET REQUIRED filesystem system timer)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Found Boost @ ${Boost_INCLUDE_DIR}")
ELSE(Boost_FOUND)
  MESSAGE(SEND_ERROR "Cannot find boost")
ENDIF(Boost_FOUND)

# CommonFile 
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/CommonFile/)

SET(ddMultigrid_SRC
    DHMTraits.h
    DHMMesh.h
    DHMMesh.cpp
    DHMField.h
    DHMField.cpp
    DHMUtil.h
    DHMUtil.cpp

    DHMAdaptiveMesh.h
    DHMAdaptiveMesh.cpp
    DHMLinearSolver.h
    DHMLinearSolver.cpp
    DHMEnergy.h
    DHMEnergy.cpp
    DHMCoarser.h
    DHMCoarser.cpp
    DHMCoarserConformal.h
    DHMCoarserConformal.cpp

    GaussLegendre.h
    GaussLegendre.cpp
#    main.cpp
)
#ADD_EXECUTABLE( ddMultigridEXEC ${ddMultigrid_SRC} )
#TARGET_LINK_LIBRARIES(ddMultigridEXEC CommonFile ${Boost_LIBRARIES})

ADD_EXECUTABLE( modal_analysis ${ddMultigrid_SRC} modal_analysis.cpp)
TARGET_LINK_LIBRARIES(modal_analysis CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( modal_analysis_time ${ddMultigrid_SRC} modal_analysis_time.cpp)
TARGET_LINK_LIBRARIES(modal_analysis_time CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( stoke ${ddMultigrid_SRC} stoke.cpp)
TARGET_LINK_LIBRARIES(stoke CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( poisson ${ddMultigrid_SRC} poisson.cpp)
TARGET_LINK_LIBRARIES(poisson CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( poisson_b ${ddMultigrid_SRC} poisson_b.cpp)
TARGET_LINK_LIBRARIES(poisson_b CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( modal_analysis_m ${ddMultigrid_SRC} modal_analysis_m.cpp)
TARGET_LINK_LIBRARIES(modal_analysis_m CommonFile ${Boost_LIBRARIES})
ADD_EXECUTABLE( stoke_m ${ddMultigrid_SRC} stoke_m.cpp)
TARGET_LINK_LIBRARIES(stoke_m CommonFile ${Boost_LIBRARIES})
