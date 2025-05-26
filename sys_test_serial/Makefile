# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2

# Executable name
EXEC = Monte_carlo_serial

# Source files
SRCS = ../serial_src/main.cpp ../serial_src/initial.cpp ../serial_src/input.cpp \
       ../serial_src/logging.cpp ../serial_src/potential.cpp ../serial_src/MC.cpp \
       ../serial_src/thermodynamic_calculator.cpp ../serial_src/rng.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Header files
HEADERS = ../serial_src/*.h

# Default target
all: $(EXEC)

# Rule to build the main executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to compile each .cpp file
../serial_src/%.o: ../serial_src/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Run simulation on all generated input files
run-inputs: $(EXEC)
	@for f in input_*.txt; do \\
		echo "Running $$f"; \\
		../$(EXEC) $$f; \\
	done

# Clean build artifacts
clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean run-inputs
