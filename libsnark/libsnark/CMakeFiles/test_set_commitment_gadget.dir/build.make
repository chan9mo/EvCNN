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
include libsnark/CMakeFiles/test_set_commitment_gadget.dir/depend.make

# Include the progress variables for this target.
include libsnark/CMakeFiles/test_set_commitment_gadget.dir/progress.make

# Include the compile flags for this target's objects.
include libsnark/CMakeFiles/test_set_commitment_gadget.dir/flags.make

libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o: libsnark/CMakeFiles/test_set_commitment_gadget.dir/flags.make
libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o: gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o -c /home/ubuntu/VCNN/libsnark/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp

libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/libsnark/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp > CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.i

libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/libsnark/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp -o CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.s

# Object files for target test_set_commitment_gadget
test_set_commitment_gadget_OBJECTS = \
"CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o"

# External object files for target test_set_commitment_gadget
test_set_commitment_gadget_EXTERNAL_OBJECTS =

libsnark/test_set_commitment_gadget: libsnark/CMakeFiles/test_set_commitment_gadget.dir/gadgetlib1/gadgets/set_commitment/tests/test_set_commitment_gadget.cpp.o
libsnark/test_set_commitment_gadget: libsnark/CMakeFiles/test_set_commitment_gadget.dir/build.make
libsnark/test_set_commitment_gadget: libsnark/libsnark.a
libsnark/test_set_commitment_gadget: depends/libff/libff/libff.a
libsnark/test_set_commitment_gadget: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/test_set_commitment_gadget: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/test_set_commitment_gadget: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libsnark/test_set_commitment_gadget: depends/libzm.a
libsnark/test_set_commitment_gadget: libsnark/CMakeFiles/test_set_commitment_gadget.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable test_set_commitment_gadget"
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_set_commitment_gadget.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsnark/CMakeFiles/test_set_commitment_gadget.dir/build: libsnark/test_set_commitment_gadget

.PHONY : libsnark/CMakeFiles/test_set_commitment_gadget.dir/build

libsnark/CMakeFiles/test_set_commitment_gadget.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/test_set_commitment_gadget.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/test_set_commitment_gadget.dir/clean

libsnark/CMakeFiles/test_set_commitment_gadget.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/test_set_commitment_gadget.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/test_set_commitment_gadget.dir/depend

