# Exact power-of-two implementation refinements

This note covers the 64-bit Q64 Falcon implementation at degrees 512/1024.
The pinned imported sources remain unchanged. `generate_pow2_backend.py` and
`generate_protocol_sampler.py` generate the adapters with asserted source
matches. The public API still performs all decoding, key binding, admission
and NTRU checks. Prepared keys and signatures are unchanged for identical
inputs. The 512 Q32 keygen and 1024 Q80 completion policies are retained.

## Trusted integer-basis preparation

`xp_pow2_decode` decodes f/g/F, recovers G, checks the bounded coefficients and
computes the unscaled public NTT key. The API verifies the stored public-key
hash before accepting the signer. `xp_pow2_prepare_core` expands that same
basis and calls `study_prepare_basis_work`. Both the trusted and general import
paths call the original basis/tree certificate and exact response-prime NTRU
check. Neither bypasses admission merely because the key was decoded.

The general importer recovers its integer basis from four inverse FFTs, then
checks four forward FFT replays. The trusted path already has the exact
integers from which its FFT basis was constructed. To show that this handoff
preserves the old check for all such inputs, let u=2^-64. The directed interval
check in `tools/check_pow2_prepare_bounds.py` verifies all 2046 nontrivial root
components within 2u of their exact roots; represented component magnitudes
are at most one. Each real Q64 product rounds with error at most u/2.

The three-real-product complex formula therefore has at most 2u rounding
error in a component. A forward butterfly consists of a root product and a
sum/difference; an inverse butterfly first sums/differences, then multiplies.
Starting with integer components of magnitude at most 127, after j stages
both exact and represented components are below

    M_j = 128 * 4^j.

Indeed represented magnitude is bounded by
`127*4^j + (4u/3)*(4^j-1) < M_j`, allowing up to 4u rounding per stage.
Writing E_j for component error against the exact-root graph, the following
recurrence covers both butterfly directions, their root errors and the
three-product rounding graph:

    E_0 = 0,
    E_(j+1) <= 4 E_j + (8 M_j + 8) u,
    E_j <= j * 4^j * 258 u.

The root error is charged against the represented operand magnitude; an
inverse butterfly can double that operand before multiplication. There are
2*(log2(N)-1) nontrivial stages across the forward/inverse pair. The final
inverse normalization divides by N/2 and adds at most u/2 rounding error.
Exact rational evaluation yields an absolute roundtrip error below 2^-24 for
both degrees. Intermediate sums and real products remain below 2^63 before
Q64 scaling, so this bound does not rely on overflowing signed carriers.

The exact transform roundtrip returns the input integer. With error below
half an integer, the original `fpr_rint` recovers that same integer. Repeating
the same forward operation graph then reproduces the stored basis bits.
Thus the trusted handoff removes checks that are necessarily satisfied by
its own bounded integer input, while the general expanded-key importer retains
all checks for externally supplied expanded data. Sign reversal and slot order
are the original mapping g,G,f,F into the reconstruction cache.

The preallocated core builds the complete expansion, mirrored widths, sample
order and response NTT cache. `XP_POW2_HANDOFF=0` keeps the old representation
validation in the same caller-owned scratch for a fair core-only control.
Wire conversion, public-key derivation, G recovery, allocation and final owner
cleanup are not part of that core timing.

## Unused sampled spectra

The existing `ffSampling_fft_inner_left_final` is selected by
`FALCON_IMPL_SKIP_DEAD_OUTPUT=1`. Integer reconstruction consumes captured
integers, not the final sampled spectra. At each node, the right subtree's
spectrum is retained because it determines the left center. The final left
subtree performs the same sample calls, but its final merges have no consumer.
Inductively, every sampled center/width, sampler call, accepted integer and
random-stream position is unchanged. Removing the unused merges cannot affect
integer reconstruction or its acceptance checks. Debug/reconstruction-check
builds retain the spectra through the existing guard.

The control is `FALCON_IMPL_SKIP_DEAD_OUTPUT=0`. The trace-only test hashes
all scalar centers, widths, outputs and complete PRNG state/cursor at each
return. All 384 attempts (589824 scalar calls) match across controls and all
three tested platforms. These traced attempts all succeeded; separate fault
tests cover the capped outer failure path. Trace instrumentation is not linked
into the normal library. Per-attempt context wiping uses volatile aligned
64-bit stores with byte edges and a compiler memory barrier.

## Q80 masked limb window

Only Falcon-1024 calls the Q80 F-only solver. Its `half_convert` needs six
adjacent base-2^31 digits of a sign-extended integer column, with three leading
zero digits. For short public lengths, the old full scan remains. For lengths
above four, `half_window` gathers the column once and uses fixed masked shift
stages to select the six digits at the secret offset.

The offset is clamped at the public sign-extension start. Each stage selects
between windows displaced by a power of two; its addresses and iteration
counts depend only on the public length. The resulting digits are exactly
those selected by the original six comparisons per limb. The final bit shift
and concatenation are unchanged. No quotient limit, rounding, reduction step,
Gaussian candidate condition or solver rejection is imported from Antrag.
Consequently all later solver decisions and seeded keys remain the same.

`XP_POW2_WINDOW=0` retains the old scan. The differential test covers 125638
cases across signed inputs, lengths 1..308, and limb/bit-shift boundaries.
Secret-taint tests cover the new window and secret scale on ARM64 and native
GCC/AVX2. The helper's bounded 316-word buffer is wiped over its entire used
prefix before return. The root/arithmetic division routines are unchanged.

## Memory, evidence and scope

The general importer uses 36*N bytes of temporary heap storage instead of its
32 KiB maximum-sized stack arrays. Its worker and wrapper now have small stack
frames. The trusted cache worker uses 20*N bytes in the core workspace, reused
after tree construction. Public preparation retains the decoded basis and
public residues alongside that workspace. All caller-owned secrets are wiped
when their owner releases them. The per-key cache and existing thread arena
have not grown.

The standalone package retains these refinements. `make check` exercises
their equivalence checks, and `VALIDATION.md` records the package validation.
These are implementation refinements, not a new Gaussian security assumption
or a full physical side-channel certification.
