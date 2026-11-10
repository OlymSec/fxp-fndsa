# Fixed-point FN-DSA: degrees 512 and 1024

A standalone C library implementing Gaussian-key Falcon with Q64.64 signing
arithmetic and the pinned experimental FN-DSA protocol. It provides key
generation, prepared signing, verification and message recovery at degrees
**512 and 1024 only**. The public API is `include/ansa/fxp_fndsa.h`.

The package targets 64-bit Linux and macOS with GCC or Clang. It contains
portable arithmetic, AVX2 kernels for native x86-64 builds, and NEON plus
runtime-checked SHA3 on AArch64. The library needs the C library and POSIX
threads; no SDK checkout, OpenSSL or GMP is required.

This is an experimental fixed-point implementation, not a claim of final
FIPS 206 compliance. The numerical admission checks and protocol scope are
specified in [PROTOCOL.md](PROTOCOL.md). The sampler uses lazy acceptance;
this package does not claim a complete constant-time or physical-side-channel
certification.

## Build and check

Install a C compiler with 128-bit integer support, GNU Make, Python 3 and
binutils (or the corresponding macOS command-line tools), then run:

```sh
make -j4
make -j4 check check-protocol check-interop check-protocol-faults
make example
```

`check-protocol-faults` uses the GNU/Linux linker's wrapping facility; omit
that target on macOS. The other tests run on both supported platforms.

For native instructions on the current host:

```sh
make BUILD=build-native NATIVE=1 -j4 check check-protocol check-interop
```

For sanitizers on Linux:

```sh
make sanitize
```

The default outputs are `build/libansa_fxp_fndsa.a` and
`build/libansa_fxp_fndsa.so.1` (macOS: `libansa_fxp_fndsa.1.dylib`).
`make install PREFIX=/your/prefix` installs the libraries and header.
`make benchmark-sdk` builds the public-API benchmark.

## API and formats

All API functions return zero on success and a negative status on failure,
except the size queries and `fxp_fndsa_free`. Size queries return zero for
unsupported degrees. Explicit key-generation and signing seeds are 64 bytes;
passing NULL uses the operating system's random source.

Prepared keys are immutable and can be shared across threads. Each signing
thread uses its own reusable scratch storage. Secret workspaces are wiped
when released. See [examples/basic.c](examples/basic.c) for an ordinary
sign/verify round trip; the header documents context, prehash and MRM calls.

| Degree | Secret key | Public key | Signature | MRM signature | Recoverable bytes |
|---:|---:|---:|---:|---:|---:|
| 512 | 1345 | 897 | 666 | 1251 | 811 |
| 1024 | 2369 | 1793 | 1280 | 2479 | 1633 |

Keys and ordinary signatures use the pinned Pornin reference's wire formats.
Ordinary signatures verify in both implementations. Seeded signatures differ
intentionally because this sampler has a distinct implementation identifier.
Accepted key-generation and preparation policies retain this fixed-point
implementation's numerical guards.

MRM uses the full fourteen-byte `Falcon MRM: H1` and `Falcon MRM: H2` tags.
The pinned reference mistakenly absorbs only two bytes; the interoperability
test checks MRM against a temporary reference object with those lengths
corrected. Recovered output is released only after authentication. An empty
recovered message is a successful result with length zero.

## Tests and source provenance

`make check` covers API errors, malformed inputs, concurrent replay, exact
NTRU relations, scalar/vector NTT equivalence, sampler replay, preparation
handoff, Q80 conversions, codecs and shared-library lifetime. `check-protocol`
independently reconstructs SHAKE transcripts and NTT equations and checks
committed test vectors. `check-interop` exercises both key generators and both
signing directions against the bundled reference. All keys in the test-vector
file are public, deterministic fixtures.

[VALIDATION.md](VALIDATION.md) records tests run for this package.
[POW2-OPTIMIZATIONS.md](POW2-OPTIMIZATIONS.md) describes the retained exact
arithmetic refinements. [provenance/package-source.json](provenance/package-source.json)
records the source commit and packaging changes. The original SDK is not
modified by this standalone release.

The ANSA authors' code is licensed under **BSD-2-Clause**. Third-party sources
retain their original licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Paper

[ansa-ibs-sc-full-v3.pdf](ansa-ibs-sc-full-v3.pdf), *Fixed-Point Antrag
Signatures and Their Application to Compact Identity-Based Signatures*, gives
concrete security bounds and a complete fixed-point implementation of Antrag,
an NTRU signature scheme with annular trapdoors, at degrees 512, 768 and 1024.
It proves ordinary and message-recovery unforgeability in the random-oracle
model, and uses message recovery in the PKG certificate to build compact
certificate-based identity-based signatures of 1,941, 2,960 and 3,884 bytes.
This library is the paper's fixed-point Falcon baseline at degrees 512 and
1024.
