# Optimized fixed-point Falcon study: SDK adapter plus unchanged vendor code.
CC ?= cc
AR ?= ar
OBJCOPY ?= objcopy
PYTHON ?= python3
BUILD ?= build
PREFIX ?= /usr/local
CFLAGS ?= -O3 -g
NATIVE ?= 0
# NTT public-key encoding, matching Antrag's two-transform verification.
NTT_PUBLIC_KEY ?= 1
ifneq ($(filter $(NTT_PUBLIC_KEY),0 1),$(NTT_PUBLIC_KEY))
$(error NTT_PUBLIC_KEY must be 0 or 1)
endif

TARGET := $(shell $(CC) -dumpmachine)
ARM64 := $(if $(filter aarch64% arm64%,$(TARGET)),1,0)
DARWIN := $(if $(findstring -apple-darwin,$(TARGET)),1,0)
ARM_SHA3 ?= $(ARM64)
# Single sequential stream; the fully unrolled XKCP C permutation wins on x64.
XKCP_SHAKE ?= $(if $(filter x86_64%,$(TARGET)),1,0)
ifneq ($(filter $(XKCP_SHAKE),0 1),$(XKCP_SHAKE))
$(error XKCP_SHAKE must be 0 or 1)
endif
V := vendor/falcon
NG := $(V)/ntrugen/src
DEFS := -DFALCON_NTT_PUBLIC_KEY=$(NTT_PUBLIC_KEY) -DFALCON_FPR_MODE=2 -DFALCON_AVX2=0 -DFALCON_PRNG_SHAKE=1 \
        -DFALCON_ARMV8_SHA3=$(ARM_SHA3) -DFALCON_MRM_EXTRACT_Q64=1 \
        -DFALCON_INTEGER_RECONSTRUCTION=1 -DFALCON_REC_SCHEDULE=1 -DFALCON_XKCP_SHAKE=$(XKCP_SHAKE)
COMMON := -std=gnu11 -fPIC -fvisibility=hidden -funroll-loops -fno-strict-aliasing \
          -Wall -Wextra -Wno-unused-function -D_DEFAULT_SOURCE -D_DARWIN_C_SOURCE
ifeq ($(NATIVE),1)
COMMON += $(if $(filter 1,$(ARM64)),-mcpu=native,-march=native)
endif
INCLUDES := -Isrc -Iinclude -I$(V) -I$(NG)
COMPILE = $(CC) $(CPPFLAGS) $(DEFS) $(COMMON) $(CFLAGS) $(INCLUDES)
FALCON_SRC := codec common falcon fft fpr rng shake sign sign_ref vrfy falcon_mrm rec_ntt exact_ntt
NG_SRC := ng_falcon ng_fxp ng_gauss ng_ntru ng_poly ng_zint31 half/solver half/division \
          $(if $(filter 1,$(ARM64)),ng_mp31_neon,ng_mp31)
VENDOR_OBJECTS := $(addprefix $(BUILD)/falcon/,$(addsuffix .o,$(FALCON_SRC))) \
                  $(addprefix $(BUILD)/ntrugen/,$(addsuffix .o,$(NG_SRC)))
ifeq ($(ARM_SHA3),1)
VENDOR_OBJECTS += $(BUILD)/falcon/shake_armv8.o
endif
ifeq ($(XKCP_SHAKE),1)
VENDOR_OBJECTS += $(BUILD)/xkcp_permute.o
endif
OBJECTS := $(BUILD)/port.o $(BUILD)/api.o $(BUILD)/protocol.o $(BUILD)/pow2_prepare.o $(VENDOR_OBJECTS)
API_HEADER := include/ansa/fxp_fndsa.h
API := $(shell $(PYTHON) tools/exports.py names $(API_HEADER))
STATIC := $(BUILD)/libansa_fxp_fndsa.a
ifeq ($(DARWIN),1)
SHARED_NAME := libansa_fxp_fndsa.1.dylib
SHARED_LINK_NAME := libansa_fxp_fndsa.dylib
SHARED_FLAGS := -dynamiclib -install_name @rpath/$(SHARED_NAME) -compatibility_version 1 \
                -Wl,-exported_symbols_list,$(BUILD)/exports.txt
localize = $(CC) -nostdlib -r -Wl,-exported_symbols_list,$(BUILD)/exports.txt -o $(2) $(1)
else
SHARED_NAME := libansa_fxp_fndsa.so.1
SHARED_LINK_NAME := libansa_fxp_fndsa.so
# The imported scratch pool registers thread-exit callbacks into this image.
SHARED_FLAGS := -shared -Wl,-z,nodelete -Wl,-soname,$(SHARED_NAME) \
                -Wl,--version-script=$(BUILD)/exports.map -Wl,--no-undefined
localize = $(CC) -nostdlib -r -o $(2) $(1) && $(OBJCOPY) $(API:%=--keep-global-symbol=%) $(2)
endif
SHARED := $(BUILD)/$(SHARED_NAME)
SHARED_LINK := $(BUILD)/$(SHARED_LINK_NAME)

.PHONY: all check check-symbols check-import sanitize install clean FORCE
all: $(STATIC) $(SHARED_LINK)
$(BUILD):
	mkdir -p $@
$(BUILD)/config: FORCE | $(BUILD)
	@printf '%s\n' '$(CC) $(CPPFLAGS) $(DEFS) $(COMMON) $(CFLAGS)' > $@.tmp
	@cmp -s $@.tmp $@ && rm $@.tmp || mv $@.tmp $@
$(BUILD)/api.o: src/api.c src/protocol.h $(API_HEADER) $(BUILD)/config
	$(CC) $(CPPFLAGS) $(DEFS) $(COMMON) $(CFLAGS) $(INCLUDES) -Werror -MMD -MP -c $< -o $@
$(BUILD)/falcon/%.o: $(V)/%.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@
$(BUILD)/ntrugen/%.o: $(NG)/%.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@
$(BUILD)/falcon/shake_armv8.o: $(V)/shake_armv8.S $(BUILD)/config
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(DEFS) -fPIC -march=armv8.2-a+sha3 -c $< -o $@
$(BUILD)/exports.map: $(API_HEADER) tools/exports.py | $(BUILD)
	$(PYTHON) tools/exports.py map $(API_HEADER) > $@
$(BUILD)/exports.txt: $(API_HEADER) tools/exports.py | $(BUILD)
	$(PYTHON) tools/exports.py darwin $(API_HEADER) > $@
$(STATIC): $(OBJECTS) $(BUILD)/exports.txt
	$(call localize,$(OBJECTS),$(BUILD)/fxp_fndsa.o)
	rm -f $@
	$(AR) rcs $@ $(BUILD)/fxp_fndsa.o
$(SHARED): $(OBJECTS) $(BUILD)/exports.map $(BUILD)/exports.txt Makefile
	$(CC) $(CFLAGS) $(SHARED_FLAGS) $(LDFLAGS) -o $@ $(OBJECTS) -pthread
$(SHARED_LINK): $(SHARED)
	ln -sf $(SHARED_NAME) $@
$(BUILD)/test_api: tests/test_api.c $(STATIC)
	$(CC) -Iinclude $(CPPFLAGS) -DFALCON_NTT_PUBLIC_KEY=$(NTT_PUBLIC_KEY) $(COMMON) $(CFLAGS) $< $(STATIC) -pthread -o $@
$(BUILD)/test_api_shared: tests/test_api.c $(SHARED_LINK)
	$(CC) -Iinclude $(CPPFLAGS) -DFALCON_NTT_PUBLIC_KEY=$(NTT_PUBLIC_KEY) $(COMMON) $(CFLAGS) $< $(SHARED_LINK) -pthread \
	    -Wl,-rpath,$(abspath $(BUILD)) -o $@
$(BUILD)/test_replay: tests/test_replay.c $(STATIC) $(VENDOR_OBJECTS) $(BUILD)/shake_original.o $(BUILD)/protocol.o
	$(COMPILE) $< $(STATIC) $(filter-out $(BUILD)/falcon/shake.o $(BUILD)/xkcp_permute.o,$(VENDOR_OBJECTS)) $(BUILD)/shake_original.o $(BUILD)/protocol.o -pthread -o $@
$(BUILD)/test_keygen: tests/test_falcon_keygen.c $(VENDOR_OBJECTS)
	$(COMPILE) $< $(VENDOR_OBJECTS) -pthread -o $@
$(BUILD)/test_keygen_api: tests/test_falcon_keygen_api.c $(VENDOR_OBJECTS)
	$(COMPILE) $< $(VENDOR_OBJECTS) -pthread -o $@
$(BUILD)/test_mq: tests/test_mq_vector.c $(VENDOR_OBJECTS)
	$(COMPILE) $< $(filter-out $(BUILD)/falcon/vrfy.o,$(VENDOR_OBJECTS)) -pthread -o $@
check-import:
	$(PYTHON) tools/check_import.py
check-symbols: all
	$(PYTHON) tools/check_symbols.py $(API_HEADER) $(STATIC) $(SHARED)
check: $(BUILD)/test_api $(BUILD)/test_api_shared $(BUILD)/test_replay \
       $(BUILD)/test_keygen $(BUILD)/test_keygen_api $(BUILD)/test_mq check-import check-symbols
	$(BUILD)/test_api
	$(BUILD)/test_api_shared
	$(BUILD)/test_replay
	$(BUILD)/test_keygen
	$(BUILD)/test_keygen_api
	$(BUILD)/test_mq
sanitize:
	$(MAKE) BUILD=build-sanitize ARM_SHA3=0 CFLAGS='-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all' check
install: all
	install -d $(DESTDIR)$(PREFIX)/include/ansa $(DESTDIR)$(PREFIX)/lib
	install -m 644 $(API_HEADER) $(DESTDIR)$(PREFIX)/include/ansa/
	install -m 644 $(STATIC) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHARED) $(DESTDIR)$(PREFIX)/lib/$(SHARED_NAME)
	ln -sf $(SHARED_NAME) $(DESTDIR)$(PREFIX)/lib/$(SHARED_LINK_NAME)
clean:
	rm -rf $(BUILD)

-include $(OBJECTS:.o=.d)

# Keep the imported Falcon SHAKE source byte-identical; generate only the
# selected permutation hook. Absorb, padding, extraction and cursors are unchanged.
ifeq ($(XKCP_SHAKE),1)
$(BUILD)/shake_backend.c: $(V)/shake.c tools/generate_shake_backend.py | $(BUILD)
	$(PYTHON) tools/generate_shake_backend.py --output $@
$(BUILD)/falcon/shake.o: $(BUILD)/shake_backend.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@
$(BUILD)/xkcp_permute.o: vendor/xkcp/keccak_permute.c vendor/xkcp/KeccakP-1600-64.macros vendor/xkcp/KeccakP-1600-unrolling.macros $(BUILD)/config
	$(COMPILE) -Ivendor/xkcp -MMD -MP -c $< -o $@
endif
$(BUILD)/shake_reference.o: $(V)/shake.c $(BUILD)/config
	$(CC) $(CPPFLAGS) $(filter-out -DFALCON_ARMV8_SHA3=%,$(DEFS)) -DFALCON_ARMV8_SHA3=0 -DFALCON_PREFIX=shake_reference $(COMMON) $(CFLAGS) $(INCLUDES) -c $< -o $@
$(BUILD)/test_shake: tests/test_shake.c $(BUILD)/shake_reference.o $(VENDOR_OBJECTS)
	$(COMPILE) $< $(BUILD)/shake_reference.o $(VENDOR_OBJECTS) -pthread -o $@
$(BUILD)/shake_original.o: $(V)/shake.c $(BUILD)/config
	$(CC) $(CPPFLAGS) $(filter-out -DFALCON_ARMV8_SHA3=%,$(DEFS)) -DFALCON_ARMV8_SHA3=0 $(COMMON) $(CFLAGS) $(INCLUDES) -c $< -o $@
.PHONY: check-shake
check-shake: $(BUILD)/test_shake
	$(PYTHON) tools/generate_shake_backend.py --check
	$(BUILD)/test_shake
check: check-shake

# Strictly matched adapters preserve the original imported source files.
$(BUILD)/falcon_backend.c: $(V)/falcon.c tools/generate_public_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_public_backend.py --kind falcon --ntt $(NTT_PUBLIC_KEY) --output $@
$(BUILD)/falcon/falcon.o: $(BUILD)/falcon_backend.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@
$(BUILD)/falcon_mrm_backend.c: $(V)/falcon_mrm.c tools/generate_public_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_public_backend.py --kind falcon_mrm --ntt $(NTT_PUBLIC_KEY) --output $@
$(BUILD)/falcon/falcon_mrm.o: $(BUILD)/falcon_mrm_backend.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@
.PHONY: check-public-backend
check-public-backend:
	$(PYTHON) tools/generate_public_backend.py --check
check: check-public-backend

$(BUILD)/protocol.o: src/protocol.c src/protocol.h $(API_HEADER) $(BUILD)/config
	$(COMPILE) -Werror -MMD -MP -c $< -o $@
$(BUILD)/sign_backend.c: $(V)/sign.c tools/generate_protocol_sampler.py $(BUILD)/capture_backend.h $(BUILD)/config
	$(PYTHON) tools/generate_protocol_sampler.py --output $@
$(BUILD)/falcon/sign.o: $(BUILD)/sign_backend.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@

# The protocol checks use independent Python hashing/encoding and a separately
# built fndsa reference in the selected build directory. No fndsa source changes.
.PHONY: check-protocol check-interop check-protocol-faults
check-protocol: all
	$(PYTHON) tests/check_protocol.py $(SHARED_LINK) --output $(BUILD)/protocol-vectors.json --expect tests/protocol-vectors.json
check-interop: all
	$(PYTHON) tests/check_fndsa_interop.py --build $(BUILD) --cc $(CC)
ifeq ($(DARWIN),0)
$(BUILD)/test_protocol_faults: tests/test_protocol_faults.c $(OBJECTS)
	$(COMPILE) $< $(sort $(OBJECTS)) -Wl,--wrap=xp_tree_attempt -Wl,--wrap=falcon_inner_get_seed -Wl,--wrap=antrag_port_scratch_push -Wl,--wrap=malloc -pthread -o $@
check-protocol-faults: $(BUILD)/test_protocol_faults
	$(BUILD)/test_protocol_faults
endif

$(BUILD)/test_protocol_codec: tests/test_protocol_codec.c $(BUILD)/protocol.o $(VENDOR_OBJECTS)
	$(COMPILE) $< $(BUILD)/protocol.o $(VENDOR_OBJECTS) -pthread -o $@
.PHONY: check-protocol-codec
check-protocol-codec: $(BUILD)/test_protocol_codec
	$(BUILD)/test_protocol_codec
check: check-protocol-codec

ifeq ($(DARWIN),0)
$(BUILD)/test_dlopen_threads: tests/test_dlopen_threads.c $(API_HEADER) | $(BUILD)
	$(CC) -Iinclude $(COMMON) $(CFLAGS) $< -ldl -pthread -o $@
check-lifecycle: $(BUILD)/test_dlopen_threads $(SHARED_LINK)
	$(BUILD)/test_dlopen_threads $(abspath $(SHARED_LINK))
check: check-lifecycle
endif

$(BUILD)/pow2_prepare.o: src/pow2_prepare.c src/pow2_prepare.h $(BUILD)/config
	$(COMPILE) -Werror -MMD -MP -c $< -o $@
$(BUILD)/capture_backend.h: $(V)/capture.h tools/generate_pow2_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_pow2_backend.py --kind capture --output $@
$(BUILD)/rec_backend.c: $(V)/rec_ntt.c tools/generate_pow2_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_pow2_backend.py --kind rec --output $@
$(BUILD)/falcon/rec_ntt.o: $(BUILD)/rec_backend.c src/pow2_clean.h $(BUILD)/config
	$(COMPILE) -MMD -MP -c $< -o $@
check-pow2-adapters:
	$(PYTHON) tools/generate_pow2_backend.py --check
check: check-pow2-adapters

$(BUILD)/half_fft_backend.h: $(NG)/half/fft.h src/pow2_window.h tools/generate_pow2_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_pow2_backend.py --kind half-fft --output $@
$(BUILD)/half_solver_backend.c: $(NG)/half/solver.c $(BUILD)/half_fft_backend.h tools/generate_pow2_backend.py $(BUILD)/config
	$(PYTHON) tools/generate_pow2_backend.py --kind half-solver --output $@
$(BUILD)/ntrugen/half/solver.o: $(BUILD)/half_solver_backend.c $(BUILD)/config
	@mkdir -p $(@D)
	$(COMPILE) -MMD -MP -c $< -o $@

$(BUILD)/test_pow2_prepare: tests/test_pow2_prepare.c $(BUILD)/pow2_prepare.o $(BUILD)/protocol.o $(VENDOR_OBJECTS)
	$(COMPILE) $^ -pthread -o $@
$(BUILD)/test_pow2_window: tests/test_pow2_window.c tests/reference_pow2_half_convert.inc $(BUILD)/half_fft_backend.h src/pow2_window.h
	$(COMPILE) -I$(BUILD) -ffunction-sections -fdata-sections $< $(if $(filter 1,$(DARWIN)),-Wl$(comma)-dead_strip,-Wl$(comma)--gc-sections) -o $@
# Trace-only build; not linked into the distributed library.
$(BUILD)/sign_trace_backend.o: $(BUILD)/sign_backend.c $(BUILD)/capture_backend.h
	$(COMPILE) -DXP_POW2_TRACE=1 -MMD -MP -c $< -o $@
$(BUILD)/test_pow2_trace: tests/test_pow2_trace.c $(BUILD)/sign_trace_backend.o $(BUILD)/protocol.o $(filter-out $(BUILD)/falcon/sign.o,$(VENDOR_OBJECTS))
	$(COMPILE) $^ -pthread -o $@
check-pow2-optimization: $(BUILD)/test_pow2_prepare $(BUILD)/test_pow2_window
	$(BUILD)/test_pow2_prepare
	$(BUILD)/test_pow2_window
check: check-pow2-optimization

comma := ,

$(BUILD)/port.o: src/port.c src/antrag_port.h $(BUILD)/config
	$(COMPILE) -MMD -MP -c $< -o $@
.PHONY: benchmark-sdk
benchmark-sdk: $(BUILD)/bench_sdk
$(BUILD)/bench_sdk: bench/bench_sdk.c $(STATIC)
	$(CC) -Iinclude $(COMMON) $(CFLAGS) $< $(STATIC) -pthread -o $@

.PHONY: example
example: $(BUILD)/example
	$(BUILD)/example
$(BUILD)/example: examples/basic.c $(STATIC)
	$(CC) -Iinclude $(CFLAGS) $< $(STATIC) -pthread -o $@
