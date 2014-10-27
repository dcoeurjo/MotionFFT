#David Coeurjolly <david.coeurjolly@liris.cnrs.fr>
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(MotionFFT)

cmake_minimum_required (VERSION 2.6)


FIND_PACKAGE(DGtal 0.8 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/")

SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies}
     ${Boost_LIBRAIRIES}
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

#fftw
FIND_PACKAGE(FFTW REQUIRED)
INCLUDE_DIRECTORIES(${FFTW_INCLUDES})



#------Tools subdirectories ------
OPTION(VERBOSE "Verbose" ON)
IF(VERBOSE)
    ADD_DEFINITIONS(-DVERBOSE)
ENDIF(VERBOSE)

SET(DGTAL_TOOLS_SRC
  createVol
  testFFT3D
  vol2raw
  )

FOREACH(FILE ${DGTAL_TOOLS_SRC})
  add_executable(${FILE} ${FILE})
target_link_libraries (${FILE} ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies} ${FFTW_LIBRARIES})
  install(TARGETS ${FILE}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
ENDFOREACH(FILE)


#------ SliceViewer -----
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
QT4_WRAP_UI(ui_sliceViewer_FORMS_HEADERS sliceViewer.ui)
QT4_WRAP_CPP(sliceViewer_HEADERS_MOC sliceViewer.h OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
ADD_EXECUTABLE(sliceViewer sliceViewer.cpp  ${sliceViewer_HEADERS_MOC}
${ui_sliceViewer_FORMS_HEADERS} )
target_link_libraries ( sliceViewer  ${DGTAL_LIBRARIES}  ${DGtalToolsLibDependencies})   
