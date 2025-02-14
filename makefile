CC = gcc
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
TEST_DIR = tests
FLAGS = -std=c99 -c -fPIC -g -Wall -include lib/qhull/src/libqhull_r/qhull_ra.h
LFLAGS = -L lib/qhull/lib -lqhull_r -lqhullstatic_r -lqhullstatic -lm
TEST_FLAGS = -g -Wall -I/usr/include/cmocka -include lib/qhull/src/libqhull_r/qhull_ra.h
OUTPUT_SO = $(BIN_DIR)/libgriddata.so

# Source files (excluding examples.c)
LIB_SRCS = $(SRC_DIR)/delaunator.c $(SRC_DIR)/interpolation.c
LIB_OBJS = $(LIB_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Example specific
EXAMPLES_SRC = $(SRC_DIR)/examples.c  # Changed from example.c to examples.c
EXAMPLE_BIN = $(BIN_DIR)/examples.exe

# Test files
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS = $(TEST_SRCS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)

# Colors for pretty output
GREEN = \033[32m
RED = \033[31m
NC = \033[0m

.PHONY: all clean test examples dirs

all: dirs lib examples

lib: dirs $(OUTPUT_SO)

examples: dirs $(EXAMPLE_BIN)

dirs:
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)

# Library target
$(OUTPUT_SO): $(LIB_OBJS)
	$(CC) -g $(LIB_OBJS) -shared -o $@ $(LFLAGS)

# Object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(FLAGS) $< -o $@

# Examples target (corrected)
$(EXAMPLE_BIN): $(EXAMPLES_SRC) $(LIB_OBJS)
	$(CC) -g $(LIB_OBJS) $< -o $@ $(LFLAGS)

# Test target
test: dirs $(TEST_BINS)
	@echo -e "$(GREEN)Running all tests...$(NC)"
	@EXIT_CODE=0; \
	echo; \
	for test in $(TEST_BINS) ; do \
		echo -e "$(GREEN)Running $$test$(NC)"; \
		if ./$$test; then \
			echo -e "$(GREEN)Test passed: $$test$(NC)"; \
		else \
			echo -e "$(RED)Test failed: $$test$(RED)"; \
			echo -e "$(NC)" \
			EXIT_CODE=1; \
		fi; \
		echo "------------------------"; \
		echo; \
	done; \
	exit $$EXIT_CODE

$(BUILD_DIR)/%: $(TEST_DIR)/%.c $(LIB_OBJS)
	$(CC) $< $(LIB_OBJS) $(TEST_FLAGS) -o $@ $(LFLAGS) -lcmocka

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

# Debug target
print-%:
	@echo $* = $($*)