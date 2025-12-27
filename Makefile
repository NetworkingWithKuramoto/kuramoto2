# Compiler
CC = gcc

# Compiler Flags
CFLAGS = -Wall -O2
CFLAGS_SERIAL = -Wall -O2 -lm 

# Linker Flags
LDFLAGS = -lm -fopenmp    # math lib + OpenMP

# Directories
SRC_DIR = src
UTILS_DIR = $(SRC_DIR)/utils
RANDOM_DIR = $(UTILS_DIR)/random
SOLVERS_DIR = $(UTILS_DIR)/solvers

# Source Files
SRC_FILES = main.c \
    $(RANDOM_DIR)/random_generators.c \
    $(SOLVERS_DIR)/gbs_integrator.c \
    $(SOLVERS_DIR)/rk4.c \
    $(SOLVERS_DIR)/rkf45.c \
    $(UTILS_DIR)/type_handlers.c \
    $(UTILS_DIR)/sampler.c \
    $(UTILS_DIR)/utils.c \
    $(UTILS_DIR)/normalizer.c

# Target Executable
# TARGET = main_paral
# Build Target Parallel
# $(TARGET): $(SRC_FILES)
# 	$(CC) $(CFLAGS) $(SRC_FILES) -o $(TARGET) $(LDFLAGS)

TARGET = main_serial
# Build Target Serial
$(TARGET): $(SRC_FILES)
	$(CC) $(SRC_FILES) -o $(TARGET) $(CFLAGS_SERIAL)

# Clean Target
clean:
	rm -f $(TARGET)
