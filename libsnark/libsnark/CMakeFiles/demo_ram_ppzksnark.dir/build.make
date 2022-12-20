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
include libsnark/CMakeFiles/demo_ram_ppzksnark.dir/depend.make

# Include the progress variables for this target.
include libsnark/CMakeFiles/demo_ram_ppzksnark.dir/progress.make

# Include the compile flags for this target's objects.
include libsnark/CMakeFiles/demo_ram_ppzksnark.dir/flags.make

libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o: libsnark/CMakeFiles/demo_ram_ppzksnark.dir/flags.make
libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o: zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o -c /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp

libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp > CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.i

libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp -o CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.s

# Object files for target demo_ram_ppzksnark
demo_ram_ppzksnark_OBJECTS = \
"CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o"

# External object files for target demo_ram_ppzksnark
demo_ram_ppzksnark_EXTERNAL_OBJECTS =

libsnark/demo_ram_ppzksnark: libsnark/CMakeFiles/demo_ram_ppzksnark.dir/zk_proof_systems/ppzksnark/ram_ppzksnark/examples/demo_ram_ppzksnark.cpp.o
libsnark/demo_ram_ppzksnark: libsnark/CMakeFiles/demo_ram_ppzksnark.dir/build.make
libsnark/demo_ram_ppzksnark: libsnark/libsnark.a
libsnark/demo_ram_ppzksnark: /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0
libsnark/demo_ram_ppzksnark: depends/libff/libff/libff.a
libsnark/demo_ram_ppzksnark: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/demo_ram_ppzksnark: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/demo_ram_ppzksnark: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libsnark/demo_ram_ppzksnark: depends/libzm.a
libsnark/demo_ram_ppzksnark: libsnark/CMakeFiles/demo_ram_ppzksnark.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable demo_ram_ppzksnark"
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo_ram_ppzksnark.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsnark/CMakeFiles/demo_ram_ppzksnark.dir/build: libsnark/demo_ram_ppzksnark

.PHONY : libsnark/CMakeFiles/demo_ram_ppzksnark.dir/build

libsnark/CMakeFiles/demo_ram_ppzksnark.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/demo_ram_ppzksnark.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/demo_ram_ppzksnark.dir/clean

libsnark/CMakeFiles/demo_ram_ppzksnark.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/demo_ram_ppzksnark.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/demo_ram_ppzksnark.dir/depend

