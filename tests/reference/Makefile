CC ?= cc
AR ?= ar
BUILD ?= build
V := vendor/fndsa_fx32_64
CFLAGS ?= -O2
SRC := $(filter-out $(V)/test_fndsa.c,$(wildcard $(V)/*.c))
OBJ := $(patsubst $(V)/%.c,$(BUILD)/%.o,$(SRC))
all: $(BUILD)/libansa_fndsa.a
$(BUILD)/%.o: $(V)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -std=gnu11 -fPIC -DFNDSA_NEON=0 -I$(V) -c $< -o $@
$(BUILD)/libansa_fndsa.a: $(OBJ)
	$(AR) rcs $@ $^
