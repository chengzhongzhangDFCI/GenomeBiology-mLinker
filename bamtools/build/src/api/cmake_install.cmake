# Install script for directory: /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/czlab/ricardo/dir_linker/dir_indel_genotyping/packages/bamtools")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/libbamtools.a")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/api_global.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAlgorithms.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAlignment.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAux.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamConstants.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamIndex.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamMultiReader.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamReader.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamWriter.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/IBamIODevice.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamConstants.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamHeader.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamProgram.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamProgramChain.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamReadGroup.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamReadGroupDictionary.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamSequence.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamSequenceDictionary.h")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/bamtools/api/algorithms" TYPE FILE FILES "/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/algorithms/Sort.h")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/internal/cmake_install.cmake")

endif()

