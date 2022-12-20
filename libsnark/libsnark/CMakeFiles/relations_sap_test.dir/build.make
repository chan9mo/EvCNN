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

# Include any dependencies generated for this target.
include libsnark/CMakeFiles/relations_sap_test.dir/depend.make

# Include the progress variables for this target.
include libsnark/CMakeFiles/relations_sap_test.dir/progress.make

# Include the compile flags for this target's objects.
include libsnark/CMakeFiles/relations_sap_test.dir/flags.make

libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o: libsnark/CMakeFiles/relations_sap_test.dir/flags.make
libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o: relations/arithmetic_programs/sap/tests/test_sap.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o -c /home/ubuntu/VCNN/libsnark/relations/arithmetic_programs/sap/tests/test_sap.cpp

libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/libsnark/relations/arithmetic_programs/sap/tests/test_sap.cpp > CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.i

libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/libsnark/relations/arithmetic_programs/sap/tests/test_sap.cpp -o CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.s

# Object files for target relations_sap_test
relations_sap_test_OBJECTS = \
"CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o"

# External object files for target relations_sap_test
relations_sap_test_EXTERNAL_OBJECTS =

libsnark/relations_sap_test: libsnark/CMakeFiles/relations_sap_test.dir/relations/arithmetic_programs/sap/tests/test_sap.cpp.o
libsnark/relations_sap_test: libsnark/CMakeFiles/relations_sap_test.dir/build.make
libsnark/relations_sap_test: libsnark/libsnark.a
libsnark/relations_sap_test: depends/libff/libff/libff.a
libsnark/relations_sap_test: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/relations_sap_test: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/relations_sap_test: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libsnark/relations_sap_test: depends/libzm.a
libsnark/relations_sap_test: libsnark/CMakeFiles/relations_sap_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable relations_sap_test"
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/relations_sap_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsnark/CMakeFiles/relations_sap_test.dir/build: libsnark/relations_sap_test

.PHONY : libsnark/CMakeFiles/relations_sap_test.dir/build

libsnark/CMakeFiles/relations_sap_test.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/relations_sap_test.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/relations_sap_test.dir/clean

libsnark/CMakeFiles/relations_sap_test.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/relations_sap_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/relations_sap_test.dir/depend

