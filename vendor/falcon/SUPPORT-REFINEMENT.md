# SHAKE and Decoder Refinement

Revision: 2026-09-07. This changes implementation mechanics only. Falcon's
parameters, fixed-point arithmetic, nonce/message transcript, SHAKE padding,
rejection threshold, signature bit order and acceptance rules are unchanged.
No prepared-verification context or public-key cache is introduced.

## ARMv8 SHA3 Permutation

Antrag's verifier uses OpenSSL SHAKE. Its native signing path separately
uses the IRS paired ARMv8 SHA3 implementation. This port takes the round
macro from `antrag_fixed_point/vendor/irs_fast/shake256x2_armv8.S`, whose
SHA-256 at the time of the port was:

```
35e00da4e12957e55a6cd49bf257649f86e929c7744ce51ce5a7cff084d19db3
```

That implementation derives from PQClean/OQS and carries an MIT license;
the notices are retained in `shake_armv8.S`. The round macro is unchanged.
Only the external state layout/load/store wrapper is adapted to one state.

Every EOR3, RAX1, XAR, BCAX and bitwise operation acts independently on
each 64-bit lane. Loading A[i] into the low half of v_i and storing that
half after 24 unchanged rounds therefore yields exactly the first lane
of the paired permutation. The round constants are Falcon's existing
standard constants. This is Keccak-f1600 on the original 25-word state,
not a new hash construction. Falcon's C permutation converts its temporary
complemented representation back before returning, so its externally stored
state is the same normal representation used by the assembly entry point.

The wrapper preserves the ABI-mandated low halves of v8 through v15 and
keeps the stack 16-byte aligned. It reads/writes exactly the 200-byte state;
the SHAKE context layout and dptr field are unchanged. Absorb, padding,
extract and value-copy semantics remain the existing implementation.

The default Linux/AArch64 build compiles only the assembly object with
the extra SHA3 instruction requirement. The C entry point checks Linux
HWCAP_SHA3, retaining the portable permutation when unavailable. A
thread-local cache stores only that public CPU capability, not hash state.
FALCON_ARMV8_SHA3=0 omits the assembly dependency entirely. Non-Linux or
non-AArch64 targets default to the existing portable path. Applications
that override general compiler architecture flags remain responsible for
their other instruction-set requirements.

## Exact Hash-to-Point Batching

The original routine reads 16-bit big-endian words, rejects w>=61445, and
reduces accepted words modulo q=12289 by repeated subtraction. Since
61445=5q, unsigned w%q is exactly the same remainder for every accepted
word. There is no approximate reciprocal or floating-point operation.

Let n be the number of coefficients still needed. Each batch extracts
2*min(n,68) bytes and processes all of them in original order. A batch
contains at most n proposals, hence cannot finish the output before its
last proposal. Rejections only increase the remaining work. Thus no unused
candidate is prefetched past the original stopping point. The routine
returns the same polynomial AND the identical complete SHAKE state/dptr
and subsequent stream suffix, not just the same output distribution.

The threshold, byte order, hash input, SHAKE rate and domain separator
are all unchanged. The constant-time hash-to-point routine is untouched;
it benefits from the same faster permutation underneath.

## MSB-First Run-Length Decoder

Antrag's LSB-first decoder cannot be copied verbatim into Falcon. The
optimization is its reservoir/run-counting strategy, translated to Falcon's
unchanged MSB-first wire format.

The new reservoir holds the next r unconsumed bits left-aligned in a
uint32_t, with all lower bits zero. A completed coefficient leaves 0<=r<=7.
One new byte supplies the next eight-bit sign/magnitude prefix: OR it at
bit 24-r, read the top eight bits, then shift left by eight. The remaining
bit count is still r. This preserves the invariant and reads exactly one
new byte, as the previous routine did.

For a nonzero reservoir, CLZ counts precisely the unary zero prefix before
the next one. This count is less than r; removing that prefix and its
terminator yields the same quotient and remaining bits. An all-zero
reservoir consumes r zeros and refills, exactly as repeated bit extraction.
CLZ is never called on zero. Shift counts are at most eight in the unary
step and 17 through 24 when appending a byte; none reaches 32.

The magnitude limit 2047, forbidden negative zero, truncation checks and
zero padding in the final consumed byte remain mandatory. Unconsumed bytes
beyond the returned length are still handled by the existing API wrapper.
Valid inputs return identical coefficients and consumed byte counts;
invalid inputs return zero. No format or verifier acceptance rule changes.

Sign restoration uses m-2*(m & (0-(s>>7))), with unsigned mask arithmetic.
Here s is 0 or 128 and 0<=m<=2047. The result is exactly m or -m and
all signed intermediates lie within [-4094,4094]. This avoids the random-sign
branch GCC generated for the previous conditional expression; it does not
change the signed coefficient or introduce a constant-time decoder claim.

## Evidence and Limits

`tests/test_support.c` checks the standard empty-input SHAKE KAT, 2000
fragmented input/output streams and context copies against the portable
permutation, and 5000 arbitrary complete Keccak states. Hash-to-point is
compared with an independent original-style two-byte/subtraction loop in
3000 cases, including rejection boundaries, complete state and suffix.

An independent bit-by-bit MSB decoder checks 22000 valid vectors over all
degrees 1 through 1024, every truncation of boundary vectors, random
truncations, mutations, arbitrary buffers, negative zero and padding.
Both actual ARM dispatch and a forced unsupported-CPU fallback are exercised
on Linux/AArch64. GCC, Clang and ASan/UBSan support tests pass. The full
fixed-Falcon suite, including NIST KATs, is also run.

These local equivalences and differential tests do not constitute formal
verification of the assembler/compiler or a new side-channel certification.
The existing variable-time hash-to-point and decoder remain variable-time.
