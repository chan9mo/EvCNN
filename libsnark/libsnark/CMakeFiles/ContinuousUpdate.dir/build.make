# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/ubuntu/VCNN

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/VCNN/libsnark

# Utility rule file for ContinuousUpdate.

# Include the progress variables for this target.
include libsnark/CMakeFiles/ContinuousUpdate.dir/progress.make

libsnark/CMakeFiles/ContinuousUpdate:
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/ctest -D ContinuousUpdate

ContinuousUpdate: libsnark/CMakeFiles/ContinuousUpdate
ContinuousUpdate: libsnark/CMakeFiles/ContinuousUpdate.dir/build.make

.PHONY : ContinuousUpdate

# Rule to build all files generated by this target.
libsnark/CMakeFiles/ContinuousUpdate.dir/build: ContinuousUpdate

.PHONY : libsnark/CMakeFiles/ContinuousUpdate.dir/build

libsnark/CMakeFiles/ContinuousUpdate.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/ContinuousUpdate.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/ContinuousUpdate.dir/clean

libsnark/CMakeFiles/ContinuousUpdate.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/ContinuousUpdate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/ContinuousUpdate.dir/depend

