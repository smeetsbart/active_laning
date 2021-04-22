#####################################
# mpacts                         #
# copyright 2020 DEM-Research Group #
#####################################


  project(laning)

  #This just checks your Cmake version 2.8 because of copy statement
  #in case of trouble 2.6 is ok if copy is omitted or replaced.
  cmake_minimum_required (VERSION 2.8)

  find_package( Boost COMPONENTS python REQUIRED )

  set( module_link_libs
    ArLite
    DEMeter
    ETility
    ${PYTHON_VERSION}
    ${Boost_PYTHON_LIBRARY}
  )
 include_directories( ${CMAKE_CURRENT_BINARY_DIR}/../include )
 link_directories( ${CMAKE_CURRENT_BINARY_DIR}/.. )

### look for Python and adjust configuration ########################
    find_package( PythonLibs )
    if( PYTHONLIBS_FOUND )
      message("Python found : ${PYTHONLIBS_VERSION_STRING}")
      include_directories( ${PYTHON_INCLUDE_DIRS} )
      message("Python libraries : ${PYTHON_LIBRARIES}")
      #link_directories( C:/src/Python27/libs )
    endif()


# make and copy the .cpp files
  file( GLOB PYTHONMODULEFILES  ${CMAKE_CURRENT_BINARY_DIR}/*.cpp)
  message( ${CMAKE_CURRENT_BINARY_DIR} ${PYTHONMODULEFILES} )
  foreach( pythonmodule ${PYTHONMODULEFILES})
    file( RELATIVE_PATH myPythonModuleRelative  ${CMAKE_CURRENT_BINARY_DIR} ${pythonmodule})
    GET_FILENAME_COMPONENT( myModuleSRCPath ${pythonmodule}           PATH    ) #source path (to find __init__.py)
    GET_FILENAME_COMPONENT( myModuleName    ${myPythonModuleRelative} NAME_WE ) #Filename without extenstion = name of module
    GET_FILENAME_COMPONENT( myModulePath    ${myPythonModuleRelative} PATH    ) #path so that we can put in the right place
    STRING(REGEX REPLACE "\\/" "." myModulePathDotted "${myModulePath}" )
    MESSAGE(STATUS "creating Python module: mpacts_usermodules${myModulePathDotted}.${myModuleName}")
  # check if  __init__ file is present and warn and ignore content when not present
    if( EXISTS ${myModuleSRCPath}/__init__.py )
    # finally make the module!
      add_library( ${myModuleName} MODULE ${pythonmodule} )
    # make sure the module is called Module and not libModule
      if(WIN32)
        set_target_properties( ${myModuleName} PROPERTIES PREFIX "" SUFFIX ".pyd"
          LIBRARY_OUTPUT_DIRECTORY_RELEASE ${myModulePath}
          LIBRARY_OUTPUT_DIRECTORY_DEBUG   ${myModulePath}
          )
      else(WIN32)
        set_target_properties( ${myModuleName} PROPERTIES PREFIX "" SUFFIX ".so"  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${myModulePath} )
      endif(WIN32)

    target_link_libraries( ${myModuleName} ${module_link_libs} )
    #done creating module!
    else( EXISTS ${myModuleSRCPath}/__init__.py )
      MESSAGE(STATUS "NO __init__.py found in ${myModuleSRCPath} doing nothing at all!")
    endif( EXISTS ${myModuleSRCPath}/__init__.py )
    #copy top level __init__.py
  endforeach(pythonmodule)

