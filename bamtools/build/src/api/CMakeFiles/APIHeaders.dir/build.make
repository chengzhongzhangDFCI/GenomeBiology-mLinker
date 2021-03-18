# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build

# Utility rule file for APIHeaders.

# Include the progress variables for this target.
include src/api/CMakeFiles/APIHeaders.dir/progress.make

APIHeaders: src/api/CMakeFiles/APIHeaders.dir/build.make
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/api_global.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/api_global.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAlgorithms.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamAlgorithms.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAlignment.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamAlignment.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamAux.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamAux.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamConstants.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamConstants.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamIndex.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamIndex.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamMultiReader.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamMultiReader.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamReader.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamReader.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/BamWriter.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/BamWriter.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/IBamIODevice.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/IBamIODevice.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamConstants.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamConstants.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamHeader.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamHeader.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamProgram.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamProgram.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamProgramChain.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamProgramChain.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamReadGroup.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamReadGroup.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamReadGroupDictionary.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamReadGroupDictionary.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamSequence.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamSequence.h
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && /usr/bin/cmake -E copy_if_different /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api/SamSequenceDictionary.h /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/include/api/SamSequenceDictionary.h
.PHONY : APIHeaders

# Rule to build all files generated by this target.
src/api/CMakeFiles/APIHeaders.dir/build: APIHeaders

.PHONY : src/api/CMakeFiles/APIHeaders.dir/build

src/api/CMakeFiles/APIHeaders.dir/clean:
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api && $(CMAKE_COMMAND) -P CMakeFiles/APIHeaders.dir/cmake_clean.cmake
.PHONY : src/api/CMakeFiles/APIHeaders.dir/clean

src/api/CMakeFiles/APIHeaders.dir/depend:
	cd /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/src/api /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api /czlab/ricardo/dir_linker/dir_indel_genotyping/bamtools/build/src/api/CMakeFiles/APIHeaders.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/api/CMakeFiles/APIHeaders.dir/depend
