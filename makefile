CC = gcc
BUILD_DIR = build
BIN_DIR = bin
SRC_DIR = src
TEST_DIR = tests
FLAGS = -std=c99 -c -fPIC -g -Wall -include lib/qhull/src/libqhull_r/qhull_ra.h
LFLAGS = -L lib/qhull/lib -lqhull_r -lqhullstatic_r -lqhullstatic -lm
TEST_FLAGS = -g -Wall -I/usr/include/cmocka -include lib/qhull/src/libqhull_r/qhull_ra.h
OUTPUT_SO = $(BIN_DIR)/libgriddata.so
OUTPUT = udi.exe

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Test files
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_OBJS = $(filter-out $(BUILD_DIR)/main.o, $(OBJS))
TEST_BINS = $(TEST_SRCS:$(TEST_DIR)/%.c=$(BUILD_DIR)/%)

# Colors for pretty output
GREEN = \033[32m
RED = \033[31m
NC = \033[0m

.PHONY: all clean test dirs

all: dirs $(OUTPUT) $(OUTPUT_SO)

dirs:
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)

$(OUTPUT): $(OBJS)
	$(CC) -g $(OBJS) -o $@ $(LFLAGS)

$(OUTPUT_SO): $(OBJS)
	$(CC) -g $(OBJS) -shared -o $@ $(LFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(FLAGS) $< -o $@

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
			EXIT_CODE=1; \
		fi; \
		echo "------------------------"; \
		echo; \
	done; \
	exit $$EXIT_CODE

$(BUILD_DIR)/%: $(TEST_DIR)/%.c $(TEST_OBJS)
	$(CC) $< $(TEST_OBJS) $(TEST_FLAGS) -o $@ $(LFLAGS) -lcmocka

clean:
	rm -rf $(OUTPUT) $(OBJS) $(OUTPUT_SO) $(BUILD_DIR) $(BIN_DIR)

# Debug target
print-%:
	@echo $* = $($*)