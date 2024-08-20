# Makefile for the unit tests

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Executable name for tests
TEST_EXEC = sys_tests

# Test source files
TEST_SRCS = sys_main_test.cpp pressure_NVT_LJ_test.cpp

# Source files from the src directory
SRCS = ../src/initial.cpp ../src/input.cpp ../src/logging.cpp ../src/potential.cpp ../src/simulation.cpp

# Object files
TEST_OBJS = $(TEST_SRCS:.cpp=.o)
SRC_OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = ../src/initial.h ../src/input.h ../src/logging.h ../src/potential.h ../src/simulation.h tests.h

# Default target (builds the test executable)
all: $(TEST_EXEC)

# Rule to build the test executable
$(TEST_EXEC): $(TEST_OBJS) $(SRC_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile test files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(TEST_OBJS) $(SRC_OBJS) $(TEST_EXEC)

# Phony targets
.PHONY: all clean
