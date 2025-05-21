# ----------- CONFIG -------------
CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
SRC_DIR = ../serial_src
TEST_SRCS = $(wildcard GCMC_ideal_particles_*.cpp)
EXECUTABLES = $(TEST_SRCS:.cpp=)

# ----------- COMMON SOURCES ------------
COMMON_SRCS = $(SRC_DIR)/MC.cpp \
              $(SRC_DIR)/cell_list.cpp \
              $(SRC_DIR)/potential.cpp \
              $(SRC_DIR)/thermodynamic_calculator.cpp \
              $(SRC_DIR)/initial.cpp \
              $(SRC_DIR)/logging.cpp \
              $(SRC_DIR)/rng.cpp \
              $(SRC_DIR)/input.cpp

# ----------- RULES --------------
all: $(EXECUTABLES)

%: %.cpp $(COMMON_SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

run: all
	@for exe in $(EXECUTABLES); do \
		echo "Running $$exe..."; \
		./$$exe & \
	done
	wait

clean:
	rm -f $(EXECUTABLES) *.o *.xyz *.txt

