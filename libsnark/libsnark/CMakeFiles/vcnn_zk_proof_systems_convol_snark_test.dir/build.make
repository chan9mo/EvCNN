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
include libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/depend.make

# Include the progress variables for this target.
include libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/progress.make

# Include the compile flags for this target's objects.
include libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/flags.make

libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o: libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/flags.make
libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o: zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o -c /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp

libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp > CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.i

libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/libsnark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/libsnark/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp -o CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.s

# Object files for target vcnn_zk_proof_systems_convol_snark_test
vcnn_zk_proof_systems_convol_snark_test_OBJECTS = \
"CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o"

# External object files for target vcnn_zk_proof_systems_convol_snark_test
vcnn_zk_proof_systems_convol_snark_test_EXTERNAL_OBJECTS =

libsnark/vcnn_zk_proof_systems_convol_snark_test: libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/zk_proof_systems/ppzksnark/convol_snark/tests/test_r1cs_vcnn_ppzksnark.cpp.o
libsnark/vcnn_zk_proof_systems_convol_snark_test: libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/build.make
libsnark/vcnn_zk_proof_systems_convol_snark_test: libsnark/libsnark.a
libsnark/vcnn_zk_proof_systems_convol_snark_test: depends/libff/libff/libff.a
libsnark/vcnn_zk_proof_systems_convol_snark_test: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/vcnn_zk_proof_systems_convol_snark_test: /usr/lib/x86_64-linux-gnu/libgmp.so
libsnark/vcnn_zk_proof_systems_convol_snark_test: /usr/lib/x86_64-linux-gnu/libgmpxx.so
libsnark/vcnn_zk_proof_systems_convol_snark_test: depends/libzm.a
libsnark/vcnn_zk_proof_systems_convol_snark_test: libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vcnn_zk_proof_systems_convol_snark_test"
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/build: libsnark/vcnn_zk_proof_systems_convol_snark_test

.PHONY : libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/build

libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/libsnark && $(CMAKE_COMMAND) -P CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/cmake_clean.cmake
.PHONY : libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/clean

libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/libsnark /home/ubuntu/VCNN/libsnark/libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsnark/CMakeFiles/vcnn_zk_proof_systems_convol_snark_test.dir/depend

