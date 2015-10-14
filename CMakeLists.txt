

cmake_minimum_required(VERSION 2.8)
 
if(POLICY CMP0020)
  cmake_policy(SET CMP0020 NEW)
endif()

PROJECT(SpuriousFibers)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

if(${VTK_VERSION} VERSION_GREATER "6" AND VTK_QT_VERSION VERSION_GREATER "4")
  # Instruct CMake to run moc automatically when needed.
  set(CMAKE_AUTOMOC ON)
  find_package(Qt5Widgets REQUIRED QUIET)
else()
  find_package(Qt4 REQUIRED)
  include(${QT_USE_FILE})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ./lib)

file(GLOB_RECURSE UI_FILES *.ui)
file(GLOB_RECURSE QT_WRAP *.h)
file(GLOB_RECURSE CXX_FILES *.cxx)

set ( QT_USE_QTXML                 TRUE )

if(${VTK_VERSION} VERSION_GREATER "6" AND VTK_QT_VERSION VERSION_GREATER "4")

  qt5_wrap_ui(UISrcs ${UI_FILES} )
  
  qt5_add_resources(ResourceSrcs ${AppResources})

  if (APPLE)
    add_executable(SpuriousFibers MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${QT_WRAP}  ${ResourceSrcs})
  else()
        add_executable(SpuriousFibers ${CXX_FILES} ${UISrcs} ${QT_WRAP}  ${ResourceSrcs})
  endif()

qt5_use_modules(SpuriousFibers Core Gui)
  
  target_link_libraries(SpuriousFibers ${VTK_LIBRARIES})

else()

  QT4_WRAP_UI(UISrcs ${UI_FILES})
  QT4_WRAP_CPP(MOCSrcs ${QT_WRAP})
  QT4_ADD_RESOURCES(ResourceSrcs ${AppResources})

  if (APPLE)
    add_executable(SpuriousFibers MACOSX_BUNDLE ${CXX_FILES} ${UISrcs} ${MOCSrcs} ${ResourceSrcs})
  else()
        add_executable(SpuriousFibers ${CXX_FILES} ${UISrcs} ${MOCSrcs} ${ResourceSrcs})
  endif()

  target_link_libraries(SpuriousFibers ${VTK_LIBRARIES} QVTK)

endif()




























