# Package validation

The standalone 512/1024 package was built and tested on Linux AArch64 with
GCC and Clang. Both runs passed `make check`, `check-protocol`, `check-interop`
and `check-protocol-faults`. The ASan/UBSan build passed `make check`, and
the OS-random key-generation/sign/verify example passed.

The 12 protocol vectors (two keys, three message modes, two degrees),
including ordinary and MRM signatures, match the selected SDK source vectors
byte for byte. Unsupported-degree size queries return zero, and degree-768
key generation returns BAD_ARGUMENT. The package contains no degree-768
sampler, key generator, object format or arithmetic implementation.

Validation covers exact NTRU equations, preparation/cache equivalence,
125,638 Q80 conversion cases, scalar/vector NTT equivalence, codecs,
malformed inputs, concurrent signing, RNG/allocation failures, capped retries,
shared-library lifetime and pinned-reference interoperability. All 17 public
symbols are checked for both static and shared libraries.

The optional interval root check in `tools/check_pow2_prepare_bounds.py`
requires Python mpmath. It is independent of the basic build/test dependencies.

`provenance/package-validation.json` records commands, compiler versions and
log digests. Linux x86-64 and macOS have not been validated locally for this
package. No new performance claim or complete
side-channel certification is made by this packaging validation.
