function(add_google_test target)
  include(GTestImport)
  import_gtest()
  target_link_libraries(${target} gtest_main gmock_main)
  add_test(NAME ${target} COMMAND ${target})
  set_property(TARGET ${target} PROPERTY FOLDER "Tests")
endfunction()

function(hivevo_add_alias target)
  # Generate aliases with the same namespace as the installed targets
  # so that they are specified with the same namespace independently of whether
  # they are imported or present in the same source tree
  add_library(hivevo::${target} ALIAS ${target})
endfunction()

function(copy_shared_lib_along_test test sharedLib)
  if (MSVC)
    #add_custom_command(TARGET ${test}
    #  POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${sharedLibFile} ${testDirectory}
    #  COMMENT "Copying shared lib alongside the test executable ${test}."
    #  VERBATIM)
    add_custom_command(TARGET ${test}
      POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${sharedLib}> $<TARGET_FILE_DIR:${test}>
      COMMENT "Copying shared lib alongside the test executable ${test}."
      VERBATIM)
  endif ()
endfunction()

function(hivevo_install_target)
####################### Install a target and mark it as exported target ######################
  set(options "")
  set(oneValueArgs EXPORT)
  set(multiValueArgs TARGETS)
  cmake_parse_arguments(HIVEVO_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # The directory "Debug" is prepended to the path for the Debug targets.
  install(TARGETS ${HIVEVO_INSTALL_TARGETS}
          EXPORT ${HIVEVO_INSTALL_EXPORT}
          LIBRARY DESTINATION $<$<CONFIG:Debug>:Debug/>lib
          ARCHIVE DESTINATION $<$<CONFIG:Debug>:Debug/>lib
          RUNTIME DESTINATION $<$<CONFIG:Debug>:Debug/>bin
          )
endfunction()

function(hivevo_install_plugin)
####################### Install a plugin and mark it as exported target ######################
  set(options "")
  set(oneValueArgs "")
  set(multiValueArgs TARGETS)
  cmake_parse_arguments(HIVEVO_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # The directory "Debug" is prepended to the path for the Debug targets.
  install(TARGETS ${HIVEVO_INSTALL_TARGETS} DESTINATION $<$<CONFIG:Debug>:Debug/>plugins)
endfunction()

function(hivevo_install_headers)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  set(options "")
  set(oneValueArgs DESTINATION)
  set(multiValueArgs DIRECTORY)
  cmake_parse_arguments(HIVEVO_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  install(DIRECTORY "${HIVEVO_INSTALL_DIRECTORY}"
          DESTINATION "${HIVEVO_INSTALL_DESTINATION}"
          COMPONENT dev
          FILES_MATCHING
          PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hxx"
          REGEX ".+/cmake" EXCLUDE  # could be omitted in new directory structure
          REGEX ".+/private" EXCLUDE # omit all private directories containing header files of pimpl classes
          )
endfunction()

function(hivevo_install_cmake_config target_name)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  # Copy the ${target_name}Config.cmake file from the source directory to the CMake binary dir
  # so that it can be installed later
  configure_file(cmake/${target_name}Config.cmake
                 "${CMAKE_CURRENT_BINARY_DIR}/${target_name}/${target_name}Config.cmake"
                 COPYONLY
                 )

  # Set destination for the generated CMake config files
  set(ConfigPackageLocation lib/cmake/${target_name})

  # Install ${target_name}Config.cmake into the specified directory
  install(FILES cmake/${target_name}Config.cmake
          DESTINATION ${ConfigPackageLocation}
          COMPONENT Devel
          )
endfunction()

function(hivevo_install_cmake_target_file target_name)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  # This makes the project importable from the build directory
  export(EXPORT ${target_name}Targets
         FILE ${CMAKE_CURRENT_BINARY_DIR}/${target_name}/${target_name}Targets.cmake
         NAMESPACE hivevo::)

  # This makes the project importable from the install directory
  # Put config file in per-project dir (name MUST match), can also
  # just go into <prefix>/cmake.
  install(EXPORT ${target_name}Targets
          FILE ${target_name}Targets.cmake
          DESTINATION lib/cmake/${target_name}
          NAMESPACE hivevo::)
endfunction()

function(hivevo_install_version_config_file target_name)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  # Generate basic XXXConfigVersion file in the CMake binary directory in the subdirectory Delib
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/${target_name}/${target_name}ConfigVersion.cmake"
                                   VERSION ${HIVEVO_VERSION}
                                   COMPATIBILITY AnyNewerVersion
                                   )

  # Install XXXConfigVersion.cmake into the specified directory
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${target_name}/${target_name}ConfigVersion.cmake"
          DESTINATION lib/cmake/${target_name}
          COMPONENT Devel
          )
endfunction()

function(hivevo_generate_export_header target_name)
  # Generate the export header file so that symbols can be exported in the shared libraries.
  include(GenerateExportHeader)
  generate_export_header(${target_name})
endfunction()

function(hivevo_install_export_header target_name)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  # the name of the export header is "lowercasetarget_export.h"
  string(TOLOWER "${target_name}" lowercase_target_name)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${lowercase_target_name}_export.h"
          DESTINATION include/${target_name}
          )
endfunction()

function(hivevo_install_dev_files)
  if (NOT HIVEVO_INSTALL_DEV_FILES)
    return()
  endif()

  set(options "")
  set(oneValueArgs DESTINATION)
  set(multiValueArgs FILES)
  cmake_parse_arguments(HIVEVO_INSTALL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  install(FILES ${HIVEVO_INSTALL_FILES}
          DESTINATION "${HIVEVO_INSTALL_DESTINATION}"
          COMPONENT dev
          )
endfunction()

function(hivevo_set_shared_library_version target_name)
  set_target_properties(${target_name} PROPERTIES
                        VERSION ${HIVEVO_VERSION}
                        SOVERSION ${HIVEVO_VERSION_MAJOR})
endfunction()

function(hivevo_install_static_target target_name)
  # Static libraries are installed only when the dev files are installed
  if (HIVEVO_INSTALL_DEV_FILES)
    hivevo_install_target(TARGETS ${target_name} EXPORT hivevoTargets)
  endif()
endfunction()

function(hivevo_install_shared_target target_name)
  # Shared libraries are always installed
  hivevo_install_target(TARGETS ${target_name} EXPORT hivevoTargets)
endfunction()

function(hivevo_install_target_dev_files target_name)
  if (HIVEVO_INSTALL_DEV_FILES)
    hivevo_install_headers(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION include/${target_name})
    hivevo_install_cmake_config(${target_name})
  endif()
endfunction()
