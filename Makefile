# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Executable name for tests
TEST_EXEC = unit_tests

# Source files for tests
TEST_SRCS = unit_test/test_main.cpp unit_test/test_energy.cpp unit_test/test_pbc.cpp

# Source files from the src directory
SRCS = src/initial.cpp src/input.cpp src/logging.cpp src/potential.cpp src/simulation.cpp

# Object files for tests
TEST_OBJS = $(TEST_SRCS:.cpp=.o)

# Object files from the src directory
OBJS = $(SRCS:.cpp=.o)

# Default target for unit tests
tests: $(TEST_EXEC)

# Rule to link object files for unit tests to create the test executable
$(TEST_EXEC): $(TEST_OBJS) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile .cpp files in unit_test directory
unit_test/%.o: unit_test/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile .cpp files in src directory
src/%.o: src/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule to remove object files and the test executable
clean:
	rm -f $(TEST_OBJS) $(OBJS) $(TEST_EXEC)

# Phony targets
.PHONY: tests clean
