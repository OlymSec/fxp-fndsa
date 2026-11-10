"""Independent SHAKE, encoding, NTT equation and MRM checks for the public profile."""
import argparse, ctypes as C, hashlib, json
from pathlib import Path
Q = 12289

class Message(C.Structure):
    _fields_ = [('data', C.c_void_p), ('data_len', C.c_size_t), ('context', C.c_void_p), ('context_len', C.c_size_t), ('hash_id', C.c_char_p)]

def load(path):
    lib = C.CDLL(str(path.resolve()))
    for n in ['public_key_bytes', 'secret_key_bytes', 'signature_bytes', 'mrm_signature_bytes', 'mrm_max_bytes']:
        f = getattr(lib, 'fxp_fndsa_' + n)
        f.argtypes = [C.c_uint]
        f.restype = C.c_size_t
    lib.fxp_fndsa_keygen.argtypes = [C.c_uint, C.c_void_p, C.c_void_p, C.c_size_t, C.c_void_p, C.c_size_t]
    lib.fxp_fndsa_prepare.argtypes = [C.POINTER(C.c_void_p), C.c_void_p, C.c_size_t]
    lib.fxp_fndsa_free.argtypes = [C.c_void_p]
    lib.fxp_fndsa_sign_message.argtypes = [C.c_void_p, C.POINTER(Message), C.c_void_p, C.c_void_p, C.c_size_t, C.POINTER(C.c_size_t)]
    lib.fxp_fndsa_verify_message.argtypes = [C.c_void_p, C.c_size_t, C.POINTER(Message), C.c_void_p, C.c_size_t]
    lib.fxp_fndsa_sign_mrm_message.argtypes = [C.c_void_p, C.c_void_p, C.c_size_t, C.POINTER(Message), C.c_void_p, C.c_void_p, C.c_size_t, C.POINTER(C.c_size_t)]
    lib.fxp_fndsa_verify_mrm_message.argtypes = [C.c_void_p, C.c_size_t, C.POINTER(Message), C.c_void_p, C.c_size_t, C.c_void_p, C.c_size_t, C.POINTER(C.c_size_t)]
    return lib

def h(x, n):
    return hashlib.shake_256(x).digest(n)

def coeffs(x, n):
    raw = h(x, 8 * n + 1024)
    out = []
    for i in range(0, len(raw), 2):
        v = int.from_bytes(raw[i:i + 2], 'little')
        if v < 61445:
            out.append(v % Q)
            if len(out) == n:
                return out
    raise AssertionError('reference XOF budget')

def fields(raw, n, bits, signed=False):
    v = int.from_bytes(raw, 'little')
    out = []
    for i in range(n):
        x = v >> i * bits & (1 << bits) - 1
        if signed and x >> bits - 1:
            x -= 1 << bits
        out.append(x)
    return out

def response(raw, n, bound):
    v = int.from_bytes(raw, 'little')
    pos = 0
    out = []
    for i in range(n):
        neg = v >> pos & 1
        mag = v >> pos + 1 & 127
        pos += 8
        while not v >> pos & 1:
            pos += 1
            mag += 128
            assert pos < len(raw) * 8 and mag <= bound
        pos += 1
        assert mag <= bound and (mag or not neg)
        out.append(-mag if neg else mag)
    assert v >> pos == 0
    return out

def evaluate(a, x):
    v = 0
    for c in reversed(a):
        v = (v * x + c) % Q
    return v

def roots(n):
    bits = n.bit_length() - 1
    psi = pow(7, 1024 // n, Q)
    return [pow(psi, 2 * int(f'{i:0{bits}b}'[::-1], 2) + 1, Q) for i in range(n)]

def transform(a):
    n = len(a)
    return [evaluate(a, x) for x in roots(n)]

def product(a, b, n):
    return [x * y % Q for (x, y) in zip(a, b)]

def inverse(a):
    n = len(a)
    out = [0] * n
    for (value, x) in zip(a, roots(n)):
        inv = pow(x, -1, Q)
        p = 1
        for j in range(n):
            out[j] = (out[j] + value * p) % Q
            p = p * inv % Q
    scale = pow(n, -1, Q)
    return [x * scale % Q for x in out]

def ring_product(a, b):
    n = len(a)
    p = [0] * (2 * n - 1)
    for (i, x) in enumerate(a):
        for (j, y) in enumerate(b):
            p[i + j] += x * y
    for i in range(2 * n - 2, n - 1, -1):
        p[i - n] -= p[i]
    return [x % Q for x in p[:n]]

def recovery_decode(v):
    pos = 0
    bits = 0
    packed = 0
    k = len(v)
    for (width, digits) in [(163, 12)] * (k // 12) + [(27, 2)] * (k % 12 // 2):
        x = sum((v[pos + i] * Q ** i for i in range(digits)))
        assert x < 1 << width
        packed |= x << bits
        bits += width
        pos += digits
    assert packed
    last = packed.bit_length() - 1
    assert last % 8 == 0
    return (packed - (1 << last)).to_bytes(last // 8, 'little')

def check_equation(n, hat, s0, s1, target):
    a = transform(s0)
    b = product(transform(s1), hat, n)
    c = transform(target)
    assert all(((x + y - z) % Q == 0 for (x, y, z) in zip(a, b, c)))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('library', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--expect', type=Path)
    a = p.parse_args()
    lib = load(a.library)
    vectors = []
    for n in [512, 1024]:
        hdr = 1
        limit = 840
        rb = {512: 625, 1024: 1239}[n]
        lam = {512: 19, 1024: 38}[n]
        gamma = {512: 15, 1024: 24}[n]
        k = n - lam - gamma
        for keyno in range(2):
            seed = bytes([n // 256, keyno]) + bytes(62)
            pk = C.create_string_buffer(lib.fxp_fndsa_public_key_bytes(n))
            sk = C.create_string_buffer(lib.fxp_fndsa_secret_key_bytes(n))
            key = C.c_void_p()
            assert lib.fxp_fndsa_keygen(n, seed, pk, len(pk), sk, len(sk)) == 0
            assert sk.raw[-64:] == h(pk.raw, 64)
            hat = fields(pk.raw[hdr:], n, 14)
            assert max(hat) < Q
            fg = 5 if n == 1024 else 6
            offset = hdr
            f = fields(sk.raw[offset:offset + n * fg // 8], n, fg, True)
            offset += n * fg // 8
            g = fields(sk.raw[offset:offset + n * fg // 8], n, fg, True)
            assert product(transform(f), hat, n) == transform(g)
            assert lib.fxp_fndsa_prepare(C.byref(key), sk, len(sk)) == 0
            try:
                for mode in ['raw', 'prehash', 'external']:
                    data = b'fixed point protocol' if mode == 'raw' else h(b'message', 64)
                    context = b'context' if mode != 'external' else b''
                    oid = b'' if mode == 'raw' else bytes.fromhex('0609608648016503040203') if mode == 'prehash' else b'\xfe'
                    db = C.create_string_buffer(data)
                    cb = C.create_string_buffer(context)
                    m = Message(C.cast(db, C.c_void_p), len(data), C.cast(cb, C.c_void_p), len(context), oid)
                    mu = data if mode == 'external' else h(h(pk.raw, 64) + bytes([int(mode == 'prehash'), len(context)]) + context + (oid if mode == 'prehash' else b'') + data, 64)
                    sig = C.create_string_buffer(lib.fxp_fndsa_mrm_signature_bytes(n))
                    written = C.c_size_t()
                    assert lib.fxp_fndsa_sign_message(key, C.byref(m), seed, sig, len(sig), C.byref(written)) == 0
                    ordinary = sig.raw[:written.value]
                    assert lib.fxp_fndsa_verify_message(pk, len(pk), C.byref(m), sig, written.value) == 0
                    s1 = response(ordinary[hdr + 40:], n, limit)
                    target = coeffs(ordinary[hdr:hdr + 40] + mu, n)
                    prod = inverse(product(transform(s1), hat, n))
                    s0 = [(x - y + Q // 2) % Q - Q // 2 for (x, y) in zip(target, prod)]
                    assert max(map(abs, s0)) <= limit
                    norm = sum((x * x for x in s0 + s1))
                    assert norm <= {512: 34034726, 1024: 70265242}[n]
                    embedded = bytes([0, 1, 0]) + h(b'recovery' + seed, lib.fxp_fndsa_mrm_max_bytes(n) - 3)
                    recovered = C.create_string_buffer(len(embedded))
                    rlen = C.c_size_t()
                    assert lib.fxp_fndsa_sign_mrm_message(key, embedded, len(embedded), C.byref(m), seed, sig, len(sig), C.byref(written)) == 0
                    mrm = sig.raw[:written.value]
                    assert lib.fxp_fndsa_verify_mrm_message(pk, len(pk), C.byref(m), sig, written.value, recovered, len(recovered), C.byref(rlen)) == 0
                    assert recovered.raw == embedded and rlen.value == len(embedded)
                    s0 = response(mrm[hdr:hdr + rb], n, limit)
                    s1 = response(mrm[hdr + rb:], n, limit)
                    target = [(x + y) % Q for (x, y) in zip(s0, inverse(product(transform(s1), hat, n)))]
                    mask = coeffs(b'Falcon MRM: H2' + b''.join((x.to_bytes(2, 'little') for x in target[:lam])), n - lam)
                    tail = [(x - y) % Q for (x, y) in zip(target[lam:], mask)]
                    rho = tail[:gamma]
                    rec = tail[gamma:]
                    assert recovery_decode(rec) == embedded
                    expected = coeffs(b'Falcon MRM: H1' + b''.join((x.to_bytes(2, 'little') for x in rho)) + mu + embedded, lam)
                    assert target[:lam] == expected
                    check_equation(n, hat, s0, s1, target)
                    norm = sum((x * x for x in s0 + s1))
                    assert norm <= {512: 34034726, 1024: 70265242}[n]
                    vectors.append({'degree': n, 'key': keyno, 'mode': mode, 'seed': seed.hex(), 'public_key': pk.raw.hex(), 'secret_key': sk.raw.hex(), 'data': data.hex(), 'context': context.hex(), 'oid': oid.hex(), 'ordinary': ordinary.hex(), 'embedded': embedded.hex(), 'mrm': mrm.hex()})
            finally:
                lib.fxp_fndsa_free(key)
        print(f'PASS {n}: independent NTT key equations, SHAKE transcripts, ordinary/MRM decoding and recovery')
    result = {'profile': 'FNDSA protocol / Q64 sampler / degrees 512 and 1024', 'vectors': vectors}
    if a.expect:
        assert result == json.loads(a.expect.read_text()), 'vector mismatch'
    a.output.write_text(json.dumps(result, indent=2) + '\n')
if __name__ == '__main__':
    main()
