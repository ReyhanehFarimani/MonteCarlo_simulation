# ----------- CONFIG -------------
CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
SRC_DIR = ../serial_src
TEST_SRCS1 = $(wildcard test_system_nvt_*.cpp)
EXECUTABLES1 = $(TEST_SRCS1:.cpp=)
TEST_SRCS2 = $(wildcard test_system_gcmc_*.cpp)
EXECUTABLES2 = $(TEST_SRCS2:.cpp=)
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
nvt: $(EXECUTABLES1)

%: %.cpp $(COMMON_SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

runnvt: nvt
	@for exe in $(EXECUTABLES1); do \
		echo "Running $$exe..."; \
		./$$exe & \
	done
	wait

cleannvt:
	rm -f $(EXECUTABLES1) *.o *.xyz *.txt

gc: $(EXECUTABLES2)

%: %.cpp $(COMMON_SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ 

rungc: gc
	@for exe in $(EXECUTABLES2); do \
		echo "Running $$exe..."; \
		./$$exe & \
	done
	wait

cleangc:
	rm -f $(EXECUTABLES2) *.o