## ======================================================================== ##
## Copyright 2009-2017 Intel Corporation                                    ##
##                                                                          ##
## Licensed under the Apache License, Version 2.0 (the "License");          ##
## you may not use this file except in compliance with the License.         ##
## You may obtain a copy of the License at                                  ##
##                                                                          ##
##     http://www.apache.org/licenses/LICENSE-2.0                           ##
##                                                                          ##
## Unless required by applicable law or agreed to in writing, software      ##
## distributed under the License is distributed on an "AS IS" BASIS,        ##
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. ##
## See the License for the specific language governing permissions and      ##
## limitations under the License.                                           ##
## ======================================================================== ##

PROJECT(ospGlutViewer)

cmake_minimum_required(VERSION 2.8)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

message("glut_libraries" ${GLUT_LIBRARIES})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


SET(OSPRAY_DIR ${CMAKE_CURRENT_LIST_DIR}/../ospray)
#SET(OSPRAY_DIR /usr/local)
FIND_LIBRARY(OSPRAY_LIBRARY
  NAMES libospray.so PATHS 
  ${OSPRAY_DIR}/clang
  ${OSPRAY_DIR}/gcc
  ${OSPRAY_DIR}/lib64 
  )
FIND_LIBRARY(OSPRAY_COMMON_LIBRARY
  NAMES libospray_common.so
  PATHS ${OSPRAY_DIR}/clang
  PATHS ${OSPRAY_DIR}/gcc
  PATHS ${OSPRAY_DIR}/lib64 
  )
# for ospray API
INCLUDE_DIRECTORIES(${OSPRAY_DIR}/ospray/include)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")

# ------------------------------------------------------------------
# take ospcommon from biff library
#
INCLUDE_DIRECTORIES(externals/biff/3rdParty/)
ADD_SUBDIRECTORY(externals/biff/3rdParty/ospcommon)

# ------------------------------------------------------------------
# now build biff
#
SET(BIFF_OSPCOMMON_LIBRARY ospray_common)
ADD_SUBDIRECTORY(externals/biff/biff)


# ------------------------------------------------------------------
# build "our" ptex
#
ADD_SUBDIRECTORY(externals/iw-ptex/)

# ------------------------------------------------------------------
# build "our" ptex
#
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/externals/biff/3rdParty #for ospcommon
  ${CMAKE_CURRENT_LIST_DIR}/externals/biff
  ${CMAKE_CURRENT_LIST_DIR}/externals/iw-ptex/src
)

ADD_SUBDIRECTORY(common/xml)
ADD_SUBDIRECTORY(common/miniSG)
ADD_SUBDIRECTORY(common/widgets)
ADD_SUBDIRECTORY(common/importer)
ADD_SUBDIRECTORY(common/commandline)

#INCLUDE_DIRECTORIES(
#  ${CMAKE_CURRENT_LIST_DIR}
#  ${CMAKE_CURRENT_LIST_DIR}/common
#  ${CMAKE_CURRENT_LIST_DIR}/../
#  ${CMAKE_CURRENT_LIST_DIR}/script
#)

# --------------------------------------------
# main executable
# --------------------------------------------

SET(APP_SRCS
  glutViewer.cpp
)

IF(OSPRAY_APPS_ENABLE_SCRIPTING)
  LIST(APPEND APP_SRCS
    ScriptedOSPGlutViewer.cpp
    ScriptedOSPGlutViewer.h
    GlutViewerScriptHandler.cpp
    GlutViewerScriptHandler.h
  )
ENDIF()

SET(LIBS
#  ospray
  ogv_commandline
  ogv_glut3d
#  ogv_minisg
#  ogv_importer
  ${OSPRAY_LIBRARY}
  ${OSPRAY_COMMON_LIBRARY}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARIES}
)

ADD_EXECUTABLE(ospGlutViewer ${APP_SRCS})
TARGET_LINK_LIBRARIES(ospGlutViewer ${LIBS})
