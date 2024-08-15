# Makefile for the main application

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Executable name
EXEC = simulation

# Source files
SRCS = src/main.cpp src/initial.cpp src/input.cpp src/logging.cpp src/potential.cpp src/simulation.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = src/initial.h src/input.h src/logging.h src/potential.h src/simulation.h

# Default target (builds the main executable)
all: $(EXEC)

# Rule to build the main executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files
src/%.o: src/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean
