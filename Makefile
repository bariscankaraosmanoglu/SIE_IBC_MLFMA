# GNU Make; building never executes the solver or regression suites.
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

CXX ?= g++
CPPFLAGS ?=
CXXFLAGS ?= -Wall -fexceptions -O3 -std=gnu++17
OPENMP_FLAGS ?= -fopenmp
LDFLAGS ?=
LDLIBS ?=
BUILD_DIR ?= build
TARGET ?= run_mlfma

# Optional GCC optimizations. Defaults retain the original release build.
# LTO=2 bounds link-time compiler jobs; ARCH/TUNE accept GCC CPU names or native.
LTO ?= 0
ARCH ?=
TUNE ?=
PGO ?= off
PROFILE_DIR ?= $(abspath $(BUILD_DIR)/profiles)

OPT_FLAGS :=
ifneq ($(LTO),0)
OPT_FLAGS += -flto=$(LTO)
endif
ifneq ($(strip $(ARCH)),)
OPT_FLAGS += -march=$(ARCH)
endif
ifneq ($(strip $(TUNE)),)
OPT_FLAGS += -mtune=$(TUNE)
endif
ifeq ($(PGO),generate)
# OpenMP workers must update profiling counters safely.
OPT_FLAGS += -fprofile-generate="$(abspath $(PROFILE_DIR))" -fprofile-update=atomic
else ifeq ($(PGO),use)
# Missing/mismatched profiles must not silently produce an untrained build.
OPT_FLAGS += -fprofile-use="$(abspath $(PROFILE_DIR))" -Werror=missing-profile -Werror=coverage-mismatch
else ifneq ($(PGO),off)
$(error PGO must be off, generate, or use)
endif

SOURCES := $(sort $(wildcard src/*.cpp))
OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPS := $(OBJECTS:.o=.d)
CONFIG := $(BUILD_DIR)/build-config.txt
# New or updated training data must invalidate a previously optimized build.
# Absolute object paths can produce nested profile directories instead of flat
# mangled names. Include both layouts without executing a shell directory walk.
profile_files = $(wildcard $(1)*.gcda) $(foreach directory,$(wildcard $(1)*/),$(call profile_files,$(directory)))
PROFILE_DATA :=
ifeq ($(PGO),use)
PROFILE_DATA := $(sort $(call profile_files,$(abspath $(PROFILE_DIR))/))
endif

# Quote settings as shell data, including embedded single quotes. Compare before
# writing so unchanged invocations retain object and executable timestamps.
shell_quote = '$(subst ','"'"',$(1))'
SETTINGS = $(foreach name,CXX CPPFLAGS CXXFLAGS OPENMP_FLAGS LDFLAGS LDLIBS LTO ARCH TUNE PGO OPT_FLAGS SOURCES PROFILE_DATA,$(call shell_quote,$(name)=$($(name))))

.PHONY: all clean help FORCE
all: $(TARGET)

$(TARGET): $(OBJECTS) $(CONFIG) Makefile
	@mkdir -p "$(@D)"
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(OPT_FLAGS) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o "$@"

$(BUILD_DIR)/%.o: src/%.cpp $(CONFIG) Makefile $(PROFILE_DATA)
	@mkdir -p "$(@D)"
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPENMP_FLAGS) $(OPT_FLAGS) -MMD -MP -MF "$(@:.o=.d)" -MT "$@" -c "$<" -o "$@"

$(CONFIG): FORCE
	@mkdir -p "$(@D)"
	@printf '%s\n' $(SETTINGS) | cmp -s - "$@" || printf '%s\n' $(SETTINGS) > "$@"

FORCE:

# Explicit opt-in only; input files, simulation outputs, and logs are retained.
clean:
	$(RM) $(OBJECTS) $(DEPS) "$(CONFIG)" "$(TARGET)"

help:
	@printf '%s\n' \
	  'make [-jN]       Build run_mlfma (no simulations or tests).' \
	  'make clean       Remove selected objects, dependencies, and executable.' \
	  'make help        Show build targets.' \
	  'Overrides: CXX CPPFLAGS CXXFLAGS OPENMP_FLAGS LDFLAGS LDLIBS BUILD_DIR TARGET' \
	  'Optional GCC switches: LTO=0 (off), LTO=2 (two LTO jobs), ARCH=native, TUNE=native' \
	  'PGO=off|generate|use; PROFILE_DIR defaults to BUILD_DIR/profiles (absolute path).' \
	  'Training runs are separate; keep the same object paths and compiler/settings for PGO use.' \
	  'Example: make -j2 LTO=2 TUNE=native BUILD_DIR=build/lto TARGET=build/lto/run_mlfma' \
	  'Example: make -j2 LTO=2 ARCH=native BUILD_DIR=build/native TARGET=build/native/run_mlfma' \
	  'Example: make -j2 BUILD_DIR=build/debug TARGET=build/debug/run_mlfma CXXFLAGS="-Wall -fexceptions -O0 -g -std=gnu++17"'

-include $(DEPS)
