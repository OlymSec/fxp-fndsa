# Build script for the Falcon implementation.
#
# ==========================(LICENSE BEGIN)============================
#
# Copyright (c) 2017-2019  Falcon Project
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ===========================(LICENSE END)=============================
#
# @author   Thomas Pornin <thomas.pornin@nccgroup.com>

.POSIX:

# =====================================================================
#
# Configurable options:
#   CC       C compiler; GCC or Clang are fine; MSVC (2015+) works too.
#   CFLAGS   Compilation flags:
#             * Optimization level -O2 or higher is recommended
#            See config.h for some possible configuration macros.
#   LD       Linker; normally the same command as the compiler.
#   LDFLAGS  Linker options, not counting the extra libs.
#   LIBS     Extra libraries for linking:
#             * If using the native FPU, test_falcon and application
#               code that calls this library may need: -lm
#               (normally not needed on x86, both 32-bit and 64-bit)

CC = clang
CFLAGS += -Wall -Wextra -Wshadow -Wundef -O3 -g #-pg -fno-pie
CPPFLAGS +=
LD = clang
LDFLAGS += #-pg -no-pie
LIBS = $(PWD)/ntrugen/src/libntrugen.a  $(PWD)/keygen_ref/build/libkeygen.a#-lm

# =====================================================================

OBJ = build/codec.o build/common.o build/falcon.o build/fft.o build/fpr.o build/rng.o build/shake.o build/sign.o build/sign_ref.o build/vrfy.o

# Only the assembly object needs SHA3 instructions; shake.c dispatches using
# Linux HWCAP and retains its portable permutation for other AArch64 CPUs.
FALCON_TARGET := $(shell $(CC) -dumpmachine)
FALCON_ARMV8_SHA3 ?= $(if $(filter aarch64% arm64%,$(FALCON_TARGET)),1,0)
ifeq ($(FALCON_ARMV8_SHA3),1)
OBJ += build/shake_armv8.o
endif
SHAKE_CPPFLAGS = -DFALCON_ARMV8_SHA3=$(FALCON_ARMV8_SHA3)

all: build build/test_falcon build/speed
all: build/test_mrm

$(OBJ) build/speed.o build/test_falcon.o build/test_mrm.o build/falcon_mrm.o build/bench_mrm.o: | build

build:
	mkdir -p build

clean:
	-rm -f $(OBJ) build/shake_armv8.o build/shake_reference.o build/test_support build/test_falcon build/test_falcon.o build/speed build/speed.o build/test_mrm build/test_mrm.o build/bench_mrm build/bench_mrm.o build/falcon_mrm.o
	-make -C ntrugen/src clean
	-make -C keygen_ref clean

test_timecop: build/test_timecop
	./timecop/run_test.sh -d 150 -s timecop/timecop.supp

prep_timecop: clean
	$(eval CFLAGS = -Wall -Wextra -Wshadow -Wundef -Wno-incompatible-pointer-types-discards-qualifiers -O0 -g)
	$(eval CPPFLAGS = -DTIMECOP)

build/test_timecop: prep_timecop ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a build/test_timecop.o build/timecop.o  $(OBJ)
	$(LD) $(LDFLAGS) -o build/test_timecop build/test_timecop.o build/timecop.o $(OBJ) $(LIBS)

build/test_falcon: ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a build/test_falcon.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/test_falcon build/test_falcon.o $(OBJ) $(LIBS)

build/test_mrm: ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a build/test_mrm.o build/falcon_mrm.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/test_mrm build/test_mrm.o build/falcon_mrm.o $(OBJ) $(LIBS)

build/bench_mrm: ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a build/bench_mrm.o build/falcon_mrm.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/bench_mrm build/bench_mrm.o build/falcon_mrm.o $(OBJ) $(LIBS)

build/speed: ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a build/speed.o $(OBJ)
	$(LD) $(LDFLAGS) -o build/speed build/speed.o $(OBJ) ntrugen/src/libntrugen.a keygen_ref/build/libkeygen.a

ntrugen/src/libntrugen.a:
	make -j -C ntrugen/src libntrugen.a

keygen_ref/build/libkeygen.a:
	CFLAGS="$(CFLAGS)" make -j -C keygen_ref

build/codec.o: codec.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/codec.o codec.c

build/common.o: common.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/common.o common.c

build/falcon.o: falcon.c falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/falcon.o falcon.c

build/fft.o: fft.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/fft.o fft.c

build/fpr.o: fpr.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/fpr.o fpr.c

build/rng.o: rng.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/rng.o rng.c

# Rebuild this small object when selecting a different SHAKE configuration.
.PHONY: force_shake_config
build/shake.o: shake.c config.h inner.h fpr.h force_shake_config
	$(CC) $(CFLAGS) $(CPPFLAGS) $(SHAKE_CPPFLAGS) -c -o build/shake.o shake.c

build/shake_armv8.o: shake_armv8.S config.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -march=armv8.2-a+sha3 -c -o $@ $<

build/sign_ref.o: sign_ref.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/sign_ref.o sign_ref.c

build/sign.o: sign.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -DFALCON_MRM_EXTRACT_Q64=1 -c -o build/sign.o sign.c

build/speed.o: speed.c falcon.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/speed.o speed.c

build/test_timecop.o: test_timecop.c falcon.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/test_timecop.o test_timecop.c

build/test_falcon.o: test_falcon.c falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/test_falcon.o test_falcon.c

build/test_mrm.o: tests/test_mrm.c falcon_mrm.h falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -I. -DFALCON_MRM_EXTRACT_Q64=1 -c -o build/test_mrm.o tests/test_mrm.c

build/bench_mrm.o: tests/bench_mrm.c falcon_mrm.h falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -I. -DFALCON_MRM_EXTRACT_Q64=1 -c -o build/bench_mrm.o tests/bench_mrm.c

build/falcon_mrm.o: falcon_mrm.c falcon_mrm.h falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -I. -DFALCON_MRM_EXTRACT_Q64=1 -c -o build/falcon_mrm.o falcon_mrm.c

build/timecop.o: timecop.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/timecop.o timecop.c

build/find_issue.o: find_issue.c falcon.h config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/find_issue.o find_issue.c

build/vrfy.o: vrfy.c config.h inner.h fpr.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o build/vrfy.o vrfy.c
