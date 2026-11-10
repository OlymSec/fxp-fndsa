# FN-DSA protocol with the Q64 Falcon sampler

This document defines the public 64-bit `fxp_fndsa_*` interface. Degrees 512
and 1024 implement the wire encodings and ordinary transcript of the pinned
Pornin `fndsa` component. MRM implements its `ds-spec.txt` with the full
specified domain tags.
The reference `fndsa` source is not modified. KRM and low-degree test profiles
are outside this interface.

## Parameters and formats

Both profiles use q=12289 and X^n+1. Hashes and streams use sequential SHAKE256. The packed bit
streams described here are least-significant bit first.

| n | Key header bytes | f/g bits | F bits | Public bytes | Secret bytes | Response bytes | Ordinary bytes | MRM bytes | Coefficient bound | Squared norm bound |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 512 | 1 | 6 | 8 | 897 | 1345 | 625 | 666 | 1251 | 840 | 34034726 |
| 1024 | 1 | 5 | 8 | 1793 | 2369 | 1239 | 1280 | 2479 | 840 | 70265242 |

At 512/1024, public-key headers are 0x09/0x0a, secret-key headers 0x59/0x5a,
ordinary-signature headers 0x39/0x3a, and MRM headers 0x29/0x2a.
A parser rejects wrong headers, wrong lengths and noncanonical encodings.

A public key contains its header and n unscaled canonical NTT residues,
14 bits each. At powers of two, entry i is h(psi^(2*rev(i)+1)), with
psi=7^(1024/n) modulo q and rev the log2(n)-bit reversal.
All residues must be below q. A verifier uses one forward and one inverse
transform without a public-key cache.

A secret key contains its header, fixed-width two's-complement f,g,F fields,
and H(vk)=SHAKE256(encoded public key)[64]. The minimum signed value of each
field is rejected; thus f/g bounds are 31 at 512 and 15 at 1024, and F is
bounded by 127. The sampler backend reconstructs G and checks its numerical
preparation requirements. The wrapper recomputes and validates the stored
H(vk). The Gaussian key-generation and admission policies remain the fixed-point
study's; this is not a claim of identical accepted-key distributions.

A response coefficient stores its sign bit, seven low magnitude bits, a run
of floor(magnitude/128) zeros, and a one terminator. Bits and bytes are LSB
first, with zero padding to the fixed response length. Negative zero, missing
terminators, excess magnitude and nonzero padding are rejected. Both response
components must obey the coefficient bound, including the ordinary
untransmitted component. Norms are squared coefficient sums.

## Message representative and ordinary signatures

The public API exposes raw-message convenience functions and `_message`
variants with context, prehash OID or external mu. For a context of at most
255 bytes:

```
H(vk) = SHAKE256(vk)[64]
pure: mu = SHAKE256(H(vk) || 0x00 || byte(ctx_len) || ctx || M)[64]
prehash: mu = SHAKE256(H(vk) || 0x01 || byte(ctx_len) || ctx || OID || digest)[64]
external: mu = caller's exact 64-byte value
```

RAW is an empty OID (NULL also means RAW). A prehash OID uses DER tag 0x06
and a short-form length at most 127; as in fndsa, digest selection and its
correct length are the caller's responsibility. EXTMU is 0xfe and ignores
context. Hash-to-coefficients reads consecutive little-endian 16-bit words,
keeps x modulo q when x<61445, and stops at the requested coefficient count.

The target is `HashToCoefficients(SHAKE256(nonce || mu), n)`, with a
40-byte nonce at every degree. An ordinary signature is header || nonce ||
encoded second response. Verification reconstructs the first response as
center(target - h*second), then checks both coefficient bounds and the norm.
No extra message framing or protocol label is inserted into this transcript.

## Hedging and retries

Let kh=SHAKE256(packed f || packed g)[40], excluding the key header and F,
as in fndsa. Let ID be eight bytes `FX64S001` for ordinary or `FX64M001` for
MRM. These identify this sampler arithmetic and must change if the same
key/message/entropy could lead to different sampler outputs in a later
implementation.

```
seedbuf = SHAKE256(ID || kh || mu || entropy)[40]
MRM additionally: seedbuf = SHAKE256(seedbuf || embedded_bytes)[40]
for counter = 0 .. 26:
    tape = SHAKE256(seedbuf || byte(counter))
    ordinary: consume 40 bytes as nonce
    MRM: consume gamma coefficients as its nonce, by the 16-bit rejection map
    derive target
    sample one candidate using the continuing tape
    if sampling, magnitude, norm, or encoding rejects: continue
    return signature
return EXHAUSTED
```

NULL seed obtains 40 OS-random bytes. The explicit test interface takes 64
bytes as entropy. The protocol changes the outer stream and target on retry;
it does not reset or reseed the stream within a candidate. The 512/1024
adapter resumes the same SHAKE state in the existing Q64 sampler.

The outer cap is a defined algorithmic limit. The fixed-point samplers retain
their stated numerical and security assumptions. The wrapper hashes f/g,
as specified by the pinned protocol.

## Message recovery

| n | lambda | gamma | recovered coefficients | Maximum recovered bytes |
|---:|---:|---:|---:|---:|
| 512 | 19 | 15 | 478 | 811 |
| 1024 | 38 | 24 | 962 | 1633 |

Compute mu from the explicit message part as above. For the recovered byte
string E, append byte 0x01 and then zeros, treating all bytes LSB first. Split
into floor(k/12) blocks of 163 bits, then (k mod 12)/2 blocks of 27 bits.
Interpret each block as a little-endian integer and emit its 12 or 2 base-q
digits in increasing significance. Decoding rejects noncanonical digits,
integers outside their 163-/27-bit domains, an all-zero bitstring, or a last
one-bit that is not at a byte boundary. The padding byte is removed.

```
c0 = HashToCoefficients(SHAKE256("Falcon MRM: H1" || u16le(rho) || mu || E), lambda)
mask = HashToCoefficients(SHAKE256("Falcon MRM: H2" || u16le(c0)), n-lambda)
c = c0 || (mask + (rho || Encode(E)) modulo q)
```

Each quoted tag has fourteen bytes, with no NUL terminator. The signature
is header || encoded first response || encoded second response. Verification
reconstructs c, removes the mask, decodes E, recomputes c0 and authenticates
before writing recovered bytes. The result status is separate from the byte
length: empty recovery returns status 0 and length 0; every error has negative
status and length 0. The pinned fndsa implementation has a two-byte-tag defect
and therefore needs a corrected reference object for these MRM vectors.

## Validation and compatibility

`tests/protocol-vectors.json` records two seeded keys at each degree with
pure, prehash and external-mu ordinary/MRM signatures. The independent Python
checker verifies the NTT key equations, transcript hashes, response norms and
MRM byte recovery. Bidirectional C tests use both key generators and ordinary
and expanded-key fndsa signing. MRM tests change only the two tag lengths in a
temporary reference object; installed fndsa files remain untouched.

The current public protocol is independent of the legacy backend's
`NTT_PUBLIC_KEY` selector. The legacy coefficient and NTT encodings remain
internal controls. The raw historical benchmark is not a measurement of
this protocol. Public API major version 1 rejects the old SDK wire contract.
