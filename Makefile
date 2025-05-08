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

# --- serial‐mode unit tests -------------------------------------
CXX         := g++
CXXFLAGS    := -std=c++17 -Iserial_src -Iunit_test_serial

# Grab all serial_src sources...
ALL_SER_SRC := $(wildcard serial_src/*.cpp)
# ...but filter out the one that has main()
SER_SRC     := $(filter-out serial_src/main.cpp, $(ALL_SER_SRC))

SER_TESTS   := $(wildcard unit_test_serial/*.cpp)
SER_BIN     := unit_test_serial/run_serial_tests

.PHONY: test_unit_serial
test_serial: $(SER_SRC) $(SER_TESTS)
	$(CXX) $(CXXFLAGS) $^ -o $(SER_BIN)
	@echo "---- Running serial unit tests ----"
	$(SER_BIN)


	# Integration‐test target
INTFLAGS    := -Iintegration_test_serial
SER_INTG    := $(wildcard integration_test_serial/*.cpp)
INTG_BIN    := integration_test_serial/run_integration_tests

.PHONY: test_integration_serial
test_integration_serial: $(SER_SRC) $(SER_INTG)
	$(CXX) $(CXXFLAGS) $(INTFLAGS) $^ -o $(INTG_BIN)
	@echo "---- Running serial integration tests ----"
	$(INTG_BIN)
