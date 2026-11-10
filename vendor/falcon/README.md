# Fixed-point C prototype for Falcon - Type 1 implementation

This folder contains the implementation of the type 1 fixed-point arithmetic stack from the paper _Toward a Secure Fixed-Point Implementation of the Falcon Signature Scheme_ (De Almeida Braga, Fouque, Lachguel, Prest).
It is based on the original Falcon reference implementation, but has been **significantly modified** to incorporate the **optimizations** and **fixed-point techniques** described in the paper.

**Disclaimer:** since this is a proof of concept implementation, it is not designed to be fully portable. 
Especially, it relies in `uint128_t` for the sake of simplicity, which is not ISO C, and likely not available on embedded devices. 
We are currently working on a better alternative implementation relying only on ISO C, memory efficient, and suited for constrained devices.

For transparency and historical context, the **original README text is preserved at the end of this file**. Option descriptions in that legacy README should be considered **legacy** unless explicitly mentioned below.

## Folder organization and key modifications

- `keygen_ref/`  
  Implements **NTRUGEN** using the **reference floating-point** implementation (native doubles, and emulated floating points mode).

- `ntrugen/`  
  Implements **NTRUGEN** using Thomas Pornin’s [fixed-point ntrugen code](https://github.com/pornin/ntrugen), with the additional bound checks and modifications described in the paper.

- Signature code (`sign.c` / `sign_ref.c`)  
  The Falcon signing path is optimized as described in the paper.  
  To facilitate comparison and differential testing, the original signing implementation was duplicated into:
  - `sign_ref.c` (reference signing path)
  - `sign.c` (optimized signing path)

- FFT / fixed-point support (`fft.c`, `fpr.c`, `fpr.h`)  
  These files contain substantial changes, including FFT adjustments and constant-time fixed-point arithmetic support required by the paper.


## Selecting floating-point vs fixed-point (compile-time)

The original README refers to `FALCON_FPNATIVE` and `FALCON_FPEMU`.  
In this codebase, the arithmetic backend is selected at compile time with:

```
FALCON_FPR_MODE={0,1,2}
  0: Native floating-point (double)
  1: Emulated floating-point (legacy; not maintained, likely broken if not used with TEST_REFERENCE)
  2: Fixed-point (default)
```

Notes:
- **Mode 0 (native)** and **Mode 2 (fixed-point)** support the optimizations described in the paper.
- **Mode 1 (emulated)** was not preserved during the refactoring and is expected to be **non-functional** outside of specific cases using the reference implementation code only (e.g., for benchmarking).


## Legacy optimization flags

The legacy implementation exposed additional compile-time flags such as `FALCON_AVX2`, `FALCON_FMA`, and `FALCON_ASM_CORTEXM4`.  
After the modifications introduced in this artifact, these flags were **not re-validated** and are **likely non-functional** at the moment.


## Build and test

Default build (fixed-point)

```sh
make
```

Native floating-point

```sh
make clean && CFLAGS="-DFALCON_FPR_MODE=0" make
```

Fixed-point (explicit)

```sh
make clean && CFLAGS="-DFALCON_FPR_MODE=2" make
```

Tests (all supported modes)

```sh
./build/test_falcon
```


## ARMv8 SHAKE and Decoder Refinements

The 2026-09-07 support-function port is documented in
`SUPPORT-REFINEMENT.md`. Linux/AArch64 builds enable the adapted Antrag/IRS
ARM SHA3 permutation, with runtime HWCAP detection and portable fallback.
`FALCON_ARMV8_SHA3=0` disables the assembly dependency. SHAKE context layout,
stream bytes, signature format and all scheme parameters remain unchanged.
Hash-to-point batches exact stream reads, and the compressed decoder uses
an MSB-first bit reservoir and run counting. No OpenSSL dependency is added.

Additional differential checks:

```sh
make -f tests/support.mk check-support
```

The controlled CPU-7 hash/decoder comparison is recorded in the enclosing
`falcon-main/benchmark-falcon-support-cpu7-20260907` directory. Its decoder
comparisons use identical coefficient arrays in the two schemes' own
unchanged formats; verification includes public-key decoding for both.

## Benchmarks

The benchmark script reproduces the four configurations used in the paper:

1. **Reference + native floating-point**
2. **Optimized + native floating-point**
3. **Reference + emulated floating-point** *(legacy; may be non-functional here)*
4. **Optimized + fixed-point**

Take note, however, that it is written to be run on a linux computer. It also tries to disable turboboost, and set the performance mode to avoid CPU frequency fluctuation, and therefore requires to run commands with root privileges. 
All benchmarked binaries are pinned to core 1.
If you need to adapt it, feel free to reuse the compilation command.

Results are stored in `./bench.log`.

```sh
./bench.sh && cat ./bench.log
```


# Original README (2020-09-30)

```
FALCON IMPLEMENTATION
=====================
Version: 2020-09-30

Falcon is a post-quantum signature algorithm, submitted to NIST's
Post-Quantum Cryptography project:

   https://csrc.nist.gov/Projects/Post-Quantum-Cryptography

Falcon is based on NTRU lattices, used with a hash-and-sign structure
and a Fourier-based sampling method that allows efficient signature
generation and verification, while producing and using relatively
compact signatures and public keys. The official Falcon Web site is:

   https://falcon-sign.info/


This implementation is written in C and is configurable at compile-time
through macros which are documented in config.h; each macro is a boolean
option and can be enabled or disabled in config.h and/or as a
command-line parameter to the compiler. Several implementation strategies
are available; however, in all cases, the same API is implemented.

Main options are the following:

  - FALCON_FPNATIVE and FALCON_FPEMU

    If using FALCON_FPNATIVE, then the C 'double' type is used for all
    floating-point operations. This is the default. This requires the
    'double' type to implement IEEE-754 semantics, in particular
    rounding to the exact precision of the 'binary64' type (i.e. "53
    bits"). The Falcon implementation takes special steps to ensure
    these properties on most common architectures. When using this
    engine, the code _may_ need to call the standard library function
    sqrt() (depending on the local architecture), which may in turn
    require linking with a specific library (e.g. adding '-lm' to the
    link command on Unix-like systems).

    FALCON_FPEMU does not use the C 'double' type, but instead works
    over only 64-bit integers and embeds its own emulation of IEEE-754
    operations. This is slower but portable, since it will work on any
    machine with a C99-compliant compiler.

  - FALCON_AVX2 and FALCON_FMA

    FALCON_AVX2, when enabled, activates the use of AVX2 compiler
    intrinsics. This works only on x86 CPU that offer AVX2 opcodes.
    Use of AVX2 improves performance. FALCON_AVX2 has no effect if
    FALCON_FPEMU is used.

    FALCON_FMA further enables the use for FMA ("fused multiply-add")
    compiler intrinsics for an extra boost to performance. This
    setting is ignored unless FALCON_FPNATIVE and FALCON_AVX2 are
    both used. Occasionally (but rarely), use of FALCON_FMA will
    change the keys and/or signatures generated from a given random
    seed, impacting reproducibility of test vectors; however, this
    has no bearing on the security of normal usage.

  - FALCON_ASM_CORTEXM4

    When enabled, inline assembly routines for FP emulation and SHAKE256
    will be used. This will work only on the ARM Cortex M3, M4 and
    compatible CPU. This assembly code is constant-time on the M4, and
    about twice faster than the generic C code used by FALCON_FPEMU.


USAGE
-----

See the Makefile for compilation flags, and config.h for configurable
options. Type 'make' to compile: this will generate two binaries called
'test_falcon' and 'speed'. 'test_falcon' runs unit tests to verify that
everything computes the expected values. 'speed' runs performance
benchmarks on Falcon-256, Falcon-512 and Falcon-1024 (Falcon-256 is a
reduced version that is faster and smaller than Falcon-512, but provides
only reduced security, and not part of the "official" Falcon).

Applications that want to use Falcon normally work on the external API,
which is documented in the "falcon.h" file. This is the only file that
an external application needs to use.

For research purposes, the inner API is documented in "inner.h". This
API gives access to many internal functions that perform some elementary
operations used in Falcon. That API also has some non-obvious
requirements, such as alignment on temporary buffers, or the need to
adjust FPU precision on 32-bit x86 systems.


LICENSE
-------

This code is provided under the MIT license:

==========================(LICENSE BEGIN)============================
Copyright (c) 2017-2020  Falcon Project

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
===========================(LICENSE END)=============================
```
