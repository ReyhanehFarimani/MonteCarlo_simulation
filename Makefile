# ==============================
# Makefile (full, object-based)
# ==============================

# -------- Main application (serial_src) --------
APP_CXX      := mpic++
APP_CXXFLAGS := -std=c++11 -Wall -Wextra -O2

EXEC := Monte_carlo_serial.o

APP_SRCS := \
  serial_src/main.cpp \
  serial_src/rng.cpp \
  serial_src/initial.cpp \
  serial_src/input.cpp \
  serial_src/logging.cpp \
  serial_src/potential.cpp \
  serial_src/thermodynamic_calculator.cpp \
  serial_src/MC.cpp \
  serial_src/cell_list.cpp \
  serial_src/GibbsMC.cpp

APP_OBJS := $(APP_SRCS:.cpp=.o)

APP_HDRS := \
  serial_src/rng.h \
  serial_src/initial.h \
  serial_src/input.h \
  serial_src/logging.h \
  serial_src/potential.h \
  serial_src/thermodynamic_calculator.h \
  serial_src/MC.h \
  serial_src/cell_list.h \
  serial_src/GibbsMC.h

.PHONY: all
all: $(EXEC)

$(EXEC): $(APP_OBJS)
	$(APP_CXX) $(APP_CXXFLAGS) -o $@ $^

serial_src/%.o: serial_src/%.cpp $(APP_HDRS)
	$(APP_CXX) $(APP_CXXFLAGS) -c $< -o $@


# -------- Serial unit tests (unit_test_serial) --------
SER_CXX      := g++
SER_CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -Iserial_src -Iunit_test_serial

# serial sources but EXCLUDE any file with '/main' in its name
SER_ALL_SRC := $(wildcard serial_src/*.cpp)
SER_SRC     := $(filter-out %/main.cpp %/main%.cpp,$(SER_ALL_SRC))

SER_TESTS   := $(wildcard unit_test_serial/*.cpp)

SER_OBJS    := $(SER_SRC:.cpp=.o) $(SER_TESTS:.cpp=.o)
SER_BIN     := unit_test_serial/run_serial_tests

.PHONY: test_serial
test_serial: $(SER_BIN)
	@echo "---- Running serial unit tests ----"
	$(SER_BIN)

$(SER_BIN): $(SER_OBJS)
	$(SER_CXX) $(SER_CXXFLAGS) $^ -o $@

serial_src/%.o: serial_src/%.cpp
	$(SER_CXX) $(SER_CXXFLAGS) -c $< -o $@
unit_test_serial/%.o: unit_test_serial/%.cpp
	$(SER_CXX) $(SER_CXXFLAGS) -c $< -o $@


# -------- Serial integration tests --------
INTFLAGS := -Iintegration_test_serial
SER_INTG := $(wildcard integration_test_serial/*.cpp)
INTG_OBJS := $(SER_INTG:.cpp=.o)
INTG_BIN  := integration_test_serial/run_integration_tests

.PHONY: test_integration_serial
test_integration_serial: $(INTG_BIN)
	@echo "---- Running serial integration tests ----"
	$(INTG_BIN)

$(INTG_BIN): $(INTG_OBJS) $(SER_SRC:.cpp=.o)
	$(SER_CXX) $(SER_CXXFLAGS) $(INTFLAGS) $^ -o $@

integration_test_serial/%.o: integration_test_serial/%.cpp
	$(SER_CXX) $(SER_CXXFLAGS) $(INTFLAGS) -c $< -o $@


# -------- Benchmarks --------
SERIAL_SRCS := $(wildcard serial_src/*.cpp)
SERIAL_OBJS := $(patsubst serial_src/%.cpp,serial_src/%.o,$(SERIAL_SRCS))
BENCH_OBJS  := $(filter-out serial_src/main.o,$(SERIAL_OBJS))

bench/benchmark_thermo: $(BENCH_OBJS) bench/benchmark_thermo.cpp
	$(SER_CXX) $(SER_CXXFLAGS) $^ -o $@

.PHONY: bench-test
bench-test: bench/benchmark_thermo
	@echo "Running benchmark..."
	@./bench/benchmark_thermo
	@echo "Benchmark completed."


# -------- MPI unit tests (unit_test_mpi) --------
MPI_CXX      := mpicxx
MPI_CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -Impi_src -Iunit_test_mpi
# Add another -I if your catch.hpp is elsewhere.

# Project MPI sources, EXCLUDE ANYTHING with '/main' in filename (robust)
MPI_ALL_SRC  := $(wildcard mpi_src/*.cpp)
MPI_LIB_SRC  := $(filter-out %/main.cpp %/main%.cpp,$(MPI_ALL_SRC))
MPI_LIB_OBJS := $(MPI_LIB_SRC:.cpp=.o)

# Tests: compile all test .cpp EXCEPT the Catch runner; add runner explicitly
MPI_TEST_ALL   := $(wildcard unit_test_mpi/*.cpp)
MPI_RUNNER_SRC := unit_test_mpi/main.cpp     # your Catch2 MPI runner
MPI_TEST_SRC   := $(filter-out $(MPI_RUNNER_SRC),$(MPI_TEST_ALL))

MPI_TEST_OBJS  := $(MPI_TEST_SRC:.cpp=.o)
MPI_RUNNER_OBJ := $(MPI_RUNNER_SRC:.cpp=.o)

MPI_BIN := unit_test_mpi/run_mpi_tests
NP     ?= 7

.PHONY: print_mpi_files
print_mpi_files:
	@echo "MPI_LIB_SRC = $(MPI_LIB_SRC)"
	@echo "MPI_TEST_SRC = $(MPI_TEST_SRC)"
	@echo "MPI_RUNNER   = $(MPI_RUNNER_SRC)"

.PHONY: test_mpi
test_mpi: $(MPI_BIN)
	@echo "---- Running MPI unit tests with $(NP) ranks ----"
	mpirun -np $(NP) $(MPI_BIN)

# Link binary from objects (each TU compiled once; only one main)
$(MPI_BIN): $(MPI_LIB_OBJS) $(MPI_TEST_OBJS) $(MPI_RUNNER_OBJ)
	$(MPI_CXX) $(MPI_CXXFLAGS) $^ -o $@

# Object build rules
mpi_src/%.o: mpi_src/%.cpp
	$(MPI_CXX) $(MPI_CXXFLAGS) -c $< -o $@

unit_test_mpi/%.o: unit_test_mpi/%.cpp
	$(MPI_CXX) $(MPI_CXXFLAGS) -c $< -o $@


# -------- Clean --------
.PHONY: clean clean_mpi clean_serial
clean: clean_mpi clean_serial
	$(RM) $(APP_OBJS) $(EXEC)

clean_serial:
	$(RM) $(SER_OBJS) $(SER_BIN) $(INTG_OBJS) $(INTG_BIN) bench/benchmark_thermo

clean_mpi:
	$(RM) $(MPI_LIB_OBJS) $(MPI_TEST_OBJS) $(MPI_RUNNER_OBJ) $(MPI_BIN)
