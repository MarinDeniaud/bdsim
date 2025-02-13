add_subdirectory(gdml)

option(USE_GZSTREAM "Zlib support through gzstream library - requires zlib" ON)

if (USE_GZSTREAM)
  message(STATUS "Looking for zlib")
  # first search for geant4 zlib
  # then general system zlib
  # then if none, disable zlib component
  if (Geant4_builtin_zlib_FOUND)
    add_definitions("-DUSE_GZSTREAM")
    if(${G4_MINOR_VERSION} GREATER 5 OR ${G4_MAJOR_VERSION} GREATER 10)
      set(ZLIB_LIBS "Geant4::G4zlib")
    else()
      set(ZLIB_LIBS "G4zlib")
    endif()
    message(STATUS "Using Geant4 built in zlib")
    add_subdirectory(gzstream)
  else()
    find_package(ZLIB)
      if(ZLIB_FOUND)
        add_definitions("-DUSE_GZSTREAM")
        set(ZLIB_LIBS ${ZLIB_LIBRARIES})
        message(STATUS "Using system zlib")
	add_subdirectory(gzstream)
      else()
        set(ZLIB_LIBS "")
        set(USE_GZSTREAM OFF)
        message(STATUS "No Geant4 included or system ZLIB found - turning off zlib support")
      endif(ZLIB_FOUND)
   endif(Geant4_builtin_zlib_FOUND)
endif(USE_GZSTREAM)
