# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2

# Executable name
EXEC = simulation

# Source files
SRCS = main.cpp initial.cpp input.cpp logging.cpp potential.cpp simulation.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = initial.h input.h logging.h potential.h simulation.h

# Default target
all: $(EXEC)

# Rule to link object files to create the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile .cpp files into .o files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule to remove object files and the executable
clean:
	rm -f $(OBJS) $(EXEC)

# Phony targets
.PHONY: all clean
