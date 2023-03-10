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
include depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/depend.make

# Include the progress variables for this target.
include depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/progress.make

# Include the compile flags for this target's objects.
include depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/flags.make

depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o: depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/flags.make
depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o: ../depends/libfqfft/tutorials/polynomial_evaluation_example.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o"
	cd /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o -c /home/ubuntu/VCNN/depends/libfqfft/tutorials/polynomial_evaluation_example.cpp

depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.i"
	cd /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/VCNN/depends/libfqfft/tutorials/polynomial_evaluation_example.cpp > CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.i

depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.s"
	cd /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/VCNN/depends/libfqfft/tutorials/polynomial_evaluation_example.cpp -o CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.s

# Object files for target polynomial_evaluation
polynomial_evaluation_OBJECTS = \
"CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o"

# External object files for target polynomial_evaluation
polynomial_evaluation_EXTERNAL_OBJECTS =

depends/libfqfft/tutorials/polynomial_evaluation: depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/polynomial_evaluation_example.cpp.o
depends/libfqfft/tutorials/polynomial_evaluation: depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/build.make
depends/libfqfft/tutorials/polynomial_evaluation: depends/libff/libff/libff.a
depends/libfqfft/tutorials/polynomial_evaluation: /usr/lib/x86_64-linux-gnu/libgmp.so
depends/libfqfft/tutorials/polynomial_evaluation: depends/libzm.a
depends/libfqfft/tutorials/polynomial_evaluation: depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/VCNN/libsnark/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable polynomial_evaluation"
	cd /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polynomial_evaluation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/build: depends/libfqfft/tutorials/polynomial_evaluation

.PHONY : depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/build

depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/clean:
	cd /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials && $(CMAKE_COMMAND) -P CMakeFiles/polynomial_evaluation.dir/cmake_clean.cmake
.PHONY : depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/clean

depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/depend:
	cd /home/ubuntu/VCNN/libsnark && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/VCNN /home/ubuntu/VCNN/depends/libfqfft/tutorials /home/ubuntu/VCNN/libsnark /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials /home/ubuntu/VCNN/libsnark/depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : depends/libfqfft/tutorials/CMakeFiles/polynomial_evaluation.dir/depend

