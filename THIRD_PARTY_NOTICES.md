# Third-party notices

The ANSA authors' code and adaptations in this distribution are licensed under
BSD-2-Clause. Third-party code retains the following licenses and notices.

| Component | Location | Origin | License |
|---|---|---|---|
| Falcon reference and fixed-point implementation | `vendor/falcon/` | Falcon Project; the fixed-point work of De Almeida Braga, Fouque, Lachguel and Prest, ePrint 2026/1531 | MIT notices in the source headers; ANSA adaptations follow the root BSD-2-Clause license |
| NTRUGEN and F-only solver | `vendor/falcon/ntrugen/` | Thomas Pornin | Unlicense in that directory; ANSA adaptations follow the root license |
| ARMv8 SHA3 permutation | `vendor/falcon/shake_armv8.S` | Adapted from PQClean/OQS; Bas Westerbaan; Hanno Becker, Vincent Hwang, Matthias J. Kannwischer, Bo-Yin Yang and Shang-Yi Yang | MIT notice in the file |
| XKCP sequential Keccak permutation | `vendor/xkcp/` | Gilles Van Assche and Ronny Van Keer; upstream commit `5585d0213203c0c5b2dfdf3b34f5da9275abc3e7` | CC0; upstream LICENSE and hashes retained |
| FN-DSA interoperability reference, tests only | `tests/reference/vendor/fndsa_fx32_64/` | Thomas Pornin, `c-fn-dsa-fx32.64`, commit `10e9c031bce4705965305bb7e5923f8cd8e286e5` | Unlicense in that directory |

The reference implementation is not linked into the distributed library.
Interoperability tests build a separate reference object; for the MRM test,
they correct the two documented domain-tag lengths in a temporary copy.

`provenance/import-sha256.json` identifies the imported optimized Falcon
snapshot. Its signing and key-generation arithmetic files are unmodified by
this packaging step. Build-time adapters remain explicit in `tools/`.
