# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -I../src -I. -Wall -Wextra

# Directories
SRC_DIR = ../src
TEST_DIR = .
BUILD_DIR = build

# Source files
TEST_SRCS = integration_test_main.cpp test_reprodicibiliy.cpp
SRC_SRCS = $(SRC_DIR)/simulation.cpp $(SRC_DIR)/input.cpp $(SRC_DIR)/logging.cpp $(SRC_DIR)/initial.cpp $(SRC_DIR)/potential.cpp 
OBJS = $(addprefix $(BUILD_DIR)/,$(TEST_SRCS:.cpp=.o)) $(SRC_SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Executable
TARGET = $(BUILD_DIR)/integration_tests

# Default target
all: $(TARGET)

# Build the test executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile test object files
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile src object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Run the tests
run: $(TARGET)
	./$(TARGET)

# Clean up build artifacts
clean:
	rm -rf $(BUILD_DIR) $(TARGET) positions1.xyz positions2.xyz simulation_data1.dat simulation_data2.dat

.PHONY: all clean run
