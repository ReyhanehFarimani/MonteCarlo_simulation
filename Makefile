# Makefile for the main application

# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Executable name
EXEC = Monte_carlo_serial

# Source files
SRCS = serial_src/main.cpp serial_src/initial.cpp serial_src/input.cpp serial_src/logging.cpp serial_src/potential.cpp serial_src/simulation.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = serial_src/initial.h serial_src/input.h serial_src/logging.h serial_src/potential.h serial_src/simulation.h

# Default target (builds the main executable)
all: $(EXEC)

# Rule to build the main executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files
serial_src/%.o: serial_src/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean
