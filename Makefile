# ==============================
# Makefile (full, object-based)
# Default build: NO PROFILING (-DNPROFILE)
# Extras: explicit *_prof targets enable profiling
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


# ======================================================
#                  M P I   S E C T I O N
# Default: NO PROFILING (adds -DNPROFILE)
# Extra explicit *_prof targets build with profiling ON
# ======================================================

# Common flags/includes
MPI_CXX      := mpicxx
MPI_INC_LIB  := -Impi_src
MPI_CXXFLAGS_BASE := -std=c++17 -O2 -Wall -Wextra $(MPI_INC_LIB)

# Profiling toggles
NPROFILE_DEF := -DNPROFILE                # disable built-in profiler
# You can also introduce: NPROFILE_HOT to disable only hot scopes if you used that macro

# ---------- MPI unit tests (unit_test_mpi) ----------
MPI_TEST_INC := -Iunit_test_mpi
MPI_CXXFLAGS := $(MPI_CXXFLAGS_BASE) $(MPI_TEST_INC) $(NPROFILE_DEF)   # DEFAULT: no profiling
MPI_CXXFLAGS_PROF := $(MPI_CXXFLAGS_BASE) $(MPI_TEST_INC)              # profiling ON

# Project MPI sources, EXCLUDE ANYTHING with '/main' in filename
MPI_ALL_SRC  := $(wildcard mpi_src/*.cpp)
MPI_LIB_SRC  := $(filter-out %/main.cpp %/main%.cpp,$(MPI_ALL_SRC))
MPI_LIB_OBJS := $(MPI_LIB_SRC:.cpp=.o)

# Tests: compile all test .cpp EXCEPT the Catch runner; add runner explicitly
MPI_TEST_ALL   := $(wildcard unit_test_mpi/*.cpp)
MPI_RUNNER_SRC := unit_test_mpi/main.cpp
MPI_TEST_SRC   := $(filter-out $(MPI_RUNNER_SRC),$(MPI_TEST_ALL))

MPI_TEST_OBJS  := $(MPI_TEST_SRC:.cpp=.o)
MPI_RUNNER_OBJ := $(MPI_RUNNER_SRC:.cpp=.o)

# Binaries
MPI_BIN_DEFAULT := unit_test_mpi/run_mpi_tests_noprof   # default
MPI_BIN_PROF    := unit_test_mpi/run_mpi_tests_prof
MPI_BIN_NOPROF  := unit_test_mpi/run_mpi_tests_noprof

# Allow overriding NP on command line: `make test_mpi NP=8`
NP ?= 4

.PHONY: print_mpi_files
print_mpi_files:
	@echo "MPI_LIB_SRC = $(MPI_LIB_SRC)"
	@echo "MPI_TEST_SRC = $(MPI_TEST_SRC)"
	@echo "MPI_RUNNER   = $(MPI_RUNNER_SRC)"

# Default test target uses NO-PROF binary
.PHONY: test_mpi
test_mpi: $(MPI_BIN_DEFAULT)
	@echo "---- Running MPI unit tests (NO PROFILING) with $(NP) ranks ----"
	mpirun -np $(NP) $(MPI_BIN_DEFAULT)

# Explicit profiling-ON test target
.PHONY: test_mpi_prof
test_mpi_prof: $(MPI_BIN_PROF)
	@echo "---- Running MPI unit tests (PROFILING ON) with $(NP) ranks ----"
	mpirun -np $(NP) $(MPI_BIN_PROF)

# Explicit profiling-OFF test target
.PHONY: test_mpi_noprof
test_mpi_noprof: $(MPI_BIN_NOPROF)
	@echo "---- Running MPI unit tests (NO PROFILING) with $(NP) ranks ----"
	mpirun -np $(NP) $(MPI_BIN_NOPROF)

# Link rules for the two unit-test executables
# NOPROF (default): objects built with -DNPROFILE
$(MPI_BIN_NOPROF): CXXFLAGS_LOCAL = $(MPI_CXXFLAGS)
$(MPI_BIN_NOPROF): $(MPI_LIB_OBJS) $(MPI_TEST_OBJS) $(MPI_RUNNER_OBJ)
	$(MPI_CXX) $(CXXFLAGS_LOCAL) $^ -o $@

# PROF: rebuild objects with profiling-enabled flags as needed
$(MPI_BIN_PROF): CXXFLAGS_LOCAL = $(MPI_CXXFLAGS_PROF)
$(MPI_BIN_PROF): $(MPI_LIB_OBJS) $(MPI_TEST_OBJS) $(MPI_RUNNER_OBJ)
	$(MPI_CXX) $(CXXFLAGS_LOCAL) $^ -o $@

# Per-directory object rules (inherit MPI_CXXFLAGS by default)
mpi_src/%.o: mpi_src/%.cpp
	$(MPI_CXX) $(MPI_CXXFLAGS) -c $< -o $@

unit_test_mpi/%.o: unit_test_mpi/%.cpp
	$(MPI_CXX) $(MPI_CXXFLAGS) -c $< -o $@


# ---------- MPI application (mpi_src) ----------
MPI_APP_CXX      := mpicxx
MPI_APP_CXXFLAGS := $(MPI_CXXFLAGS_BASE) $(NPROFILE_DEF)   # DEFAULT APP: NO PROFILING

# Path to the MPI main translation unit (override if needed)
# e.g.: make mpi_app MPI_APP_MAIN=mpi_src/main_mpi.cpp
MPI_APP_MAIN ?= mpi_src/main.cpp

# Reuse library objects from the MPI unit-test section (all mpi_src/*.cpp except main*)
MPI_APP_BIN_DEFAULT := Monte_carlo_mpi_noprof   # default binary (no profiling)
MPI_APP_BIN_PROF    := Monte_carlo_mpi_prof     # profiling ON
MPI_APP_BIN_NOPROF  := Monte_carlo_mpi_noprof   # explicit no-profiling

# Default app target (NO PROFILING)
.PHONY: mpi_app
mpi_app: $(MPI_APP_BIN_DEFAULT)

# Build default (no-profiling) app
$(MPI_APP_BIN_DEFAULT): $(MPI_LIB_OBJS) $(MPI_APP_MAIN:.cpp=.o)
	$(MPI_APP_CXX) $(MPI_APP_CXXFLAGS) $^ -o $@

# Ensure we can compile the main TU (inherits MPI_APP_CXXFLAGS which contain -DNPROFILE by default)
$(MPI_APP_MAIN:.cpp=.o): $(MPI_APP_MAIN)
	$(MPI_APP_CXX) $(MPI_APP_CXXFLAGS) -c $< -o $@

# Explicit profiling-ON app
.PHONY: mpi_app_prof
mpi_app_prof: $(MPI_APP_BIN_PROF)

$(MPI_APP_BIN_PROF): override MPI_APP_CXXFLAGS := $(MPI_CXXFLAGS_BASE)    # profiling ON
$(MPI_APP_BIN_PROF): $(MPI_LIB_OBJS) $(MPI_APP_MAIN:.cpp=.o)
	$(MPI_APP_CXX) $(MPI_APP_CXXFLAGS) $^ -o $@

# Explicit no-profiling app (same as default; kept for symmetry)
.PHONY: mpi_app_noprof
mpi_app_noprof: $(MPI_APP_BIN_NOPROF)

$(MPI_APP_BIN_NOPROF): override MPI_APP_CXXFLAGS := $(MPI_CXXFLAGS_BASE) $(NPROFILE_DEF)
$(MPI_APP_BIN_NOPROF): $(MPI_LIB_OBJS) $(MPI_APP_MAIN:.cpp=.o)
	$(MPI_APP_CXX) $(MPI_APP_CXXFLAGS) $^ -o $@


# -------- Sanity-check mini apps (parallel + serial) --------
# Parallel sanity uses the default NO-PROF flags via MPI_CXXFLAGS (= ... -DNPROFILE)

# Compilers/flags (reuse your project include dirs)
SAN_SER_CXX      := g++
SAN_SER_CXXFLAGS := -std=c++17 -O3 -Wall -Wextra -Iserial_src -Isanity_check

SAN_MPI_CXX      := mpicxx
SAN_MPI_CXXFLAGS := $(MPI_CXXFLAGS) -Isanity_check

# Binaries
SAN_PAR_BIN := sanity_check/run_parallel
SAN_SER_BIN := sanity_check/run_serial

# Mains
SAN_PAR_MAIN := sanity_check/main_parallel.cpp
SAN_SER_MAIN := sanity_check/main_serial.cpp

# Objects for mains (compiled with the correct compilers)
SAN_PAR_MAIN_OBJ := sanity_check/main_parallel.o
SAN_SER_MAIN_OBJ := sanity_check/main_serial.o

# Reuse your existing project object lists so we link against the same lib code:
SAN_SER_LIB_OBJS := $(SER_SRC:.cpp=.o)

# Build both sanity apps
.PHONY: sanity_build
sanity_build: $(SAN_PAR_BIN) $(SAN_SER_BIN)

# Parallel sanity exe: link MPI library objects + parallel main object (NO PROFILING by default)
$(SAN_PAR_BIN): $(MPI_LIB_OBJS) $(SAN_PAR_MAIN_OBJ)
	$(SAN_MPI_CXX) $(SAN_MPI_CXXFLAGS) $^ -o $@

# Serial sanity exe: link serial library objects + serial main object
$(SAN_SER_BIN): $(SAN_SER_LIB_OBJS) $(SAN_SER_MAIN_OBJ)
	$(SAN_SER_CXX) $(SAN_SER_CXXFLAGS) $^ -o $@

# Compile sanity-check mains with the correct compilers
sanity_check/main_parallel.o: sanity_check/main_parallel.cpp
	$(SAN_MPI_CXX) $(SAN_MPI_CXXFLAGS) -c $< -o $@

sanity_check/main_serial.o: sanity_check/main_serial.cpp
	$(SAN_SER_CXX) $(SAN_SER_CXXFLAGS) -c $< -o $@

# Convenience: run both apps and the Python comparator
NP ?= 4
SAN_BASE ?= sanity_check/out

.PHONY: sanity_run
sanity_run: sanity_build
	@echo "---- Running sanity parallel ($(NP) ranks) ----"
	mpirun -np $(NP) $(SAN_PAR_BIN) $(SAN_BASE)
	@echo "---- Running sanity serial ----"
	$(SAN_SER_BIN) $(SAN_BASE)
	@echo "---- Python comparison ----"
	python3 sanity_check/check_sanity.py $(SAN_BASE)


# -------- Clean --------
.PHONY: clean clean_mpi clean_serial clean_sanity
clean: clean_mpi clean_serial clean_sanity
	$(RM) $(APP_OBJS) $(EXEC)

clean_serial:
	$(RM) $(SER_OBJS) $(SER_BIN) $(INTG_OBJS) $(INTG_BIN) bench/benchmark_thermo

clean_mpi:
	$(RM) $(MPI_LIB_OBJS) $(MPI_TEST_OBJS) $(MPI_RUNNER_OBJ) \
	      $(MPI_BIN_DEFAULT) $(MPI_BIN_PROF) $(MPI_BIN_NOPROF) \
	      $(MPI_APP_BIN_DEFAULT) $(MPI_APP_BIN_PROF) $(MPI_APP_BIN_NOPROF) \
	      $(MPI_APP_MAIN:.cpp=.o)

clean_sanity:
	$(RM) $(SAN_PAR_MAIN_OBJ) $(SAN_SER_MAIN_OBJ) $(SAN_PAR_BIN) $(SAN_SER_BIN)
