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
include libsnark/CMakeFiles/new_commit_vcnn.dir/depend.make

# Include the progress variables for this target.
include libsnark/CMakeFiles/new_commit_vcnn.dir/progress.make

# Include the compile flags for this target's objects.
include libsnark/CMakeFiles/new_commit_vcnn.dir/flags.make

libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o: libsnark/CMakeFiles/new_commit_vcnn.dir/flags.make
libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o: zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o -c /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp

libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp > CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i

libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp -o CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s

# Object files for target new_commit_vcnn
new_commit_vcnn_OBJECTS = \
"CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o"

# External object files for target new_commit_vcnn
new_commit_vcnn_EXTERNAL_OBJECTS =

libsnark/new_commit_vcnn: libsnark/CMakeFiles/new_commit_vcnn.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o
libsnark/new_commit_vcnn: libsnark/CMakeFiles/new_commit_vcnn.dir/build.make
libsnark/new_commit_vcnn: libsnark/libsnark.a
libsnark/new_commit_vcnn: depends/libff/libff/libff.a
libsnark/new_commit_vcnn: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/new_commit_vcnn: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/new_commit_vcnn: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libsnark/new_commit_vcnn: depends/libzm.a
libsnark/new_commit_vcnn: libsnark/CMakeFiles/new_commit_vcnn.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable new_commit_vcnn"
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/new_commit_vcnn.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsnark/CMakeFiles/new_commit_vcnn.dir/build: libsnark/new_commit_vcnn

.PHONY : libsnark/CMakeFiles/new_commit_vcnn.dir/build

libsnark/CMakeFiles/new_commit_vcnn.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/new_commit_vcnn.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/new_commit_vcnn.dir/clean

libsnark/CMakeFiles/new_commit_vcnn.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/new_commit_vcnn.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/new_commit_vcnn.dir/depend
