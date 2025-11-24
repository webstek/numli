
# Simple Makefile to build tests into executables

CXX := g++
# Optional AVX2 build: set `AVX2=1` when invoking make (e.g. `make AVX2=1`)
AVX2 ?= 0
ifeq ($(AVX2),1)
AVX2_FLAGS := -mavx2 -mbmi2 -mfma
else
AVX2_FLAGS :=
endif

CXXFLAGS := -std=c++23 -g -Wall -Wextra -I. -Isrc -Ivendor -Itests $(AVX2_FLAGS)

SRC_DIR := tests
BUILD_DIR := build

SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
TARGETS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%,$(SOURCES))

all: $(TARGETS)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

run-%: $(BUILD_DIR)/%
	@./$(BUILD_DIR)/$*

avx2:
	@$(MAKE) --no-print-directory AVX2=1

clean:
	@rm -rf $(BUILD_DIR)/*.exe $(BUILD_DIR)/*

.PHONY: all clean run-% avx2