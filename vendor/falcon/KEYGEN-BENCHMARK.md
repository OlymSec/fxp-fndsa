# Fixed-Point Falcon Key-Generation Benchmark

Date: 2026-09-04 (Asia/Manila)

## Results

The benchmark was pinned to logical CPU 7.

| Parameter set | Time per keygen | Throughput |
| --- | ---: | ---: |
| Falcon-512 | 2.45 ms | 408.2 keygens/s |
| Falcon-1024 | 14.68 ms | 68.1 keygens/s |

## Configuration

- Architecture: AArch64, 8 logical CPUs
- Operating system: Linux 5.15.0-190-generic
- Compiler: Ubuntu Clang 14.0.0
- Arithmetic mode: fixed-point (`FALCON_FPR_MODE=2`)
- Optimization: `-O3`
- CPU affinity: logical CPU 7 only
- Timing threshold: at least 2 seconds per adaptive benchmark batch

The benchmark uses `clock()` to measure process CPU time. Each operation is
warmed up with five executions before the adaptive timed batch. CPU frequency
and governor settings were not changed for this run.

## Commands

```sh
make clean
make CFLAGS="-DFALCON_FPR_MODE=2 -Wall -Wextra -Wshadow -Wundef -O3 -g" build/speed
taskset -c 7 ./build/speed 2
```

## Raw Output

The `kg(ms)` column is the key-generation measurement recorded above. The
remaining columns are incidental measurements emitted by the stock benchmark.

```text
time threshold = 2.0000 s
kg = keygen, ek = expand private key,
st = sign (with expanded key), vv = verify
stc, vvc: like sd, st and vv, but with constant-time hash-to-point
keygen in milliseconds, other values in microseconds

degree  kg(ms)   ek(us)   st(us)  stc(us)   vv(us)  vvc(us)
 512:     2.45   104.39   138.81   147.80     9.88    16.79
1024:    14.68   225.88   300.01   321.01    21.24    34.51
```

This record contains one invocation of the stock adaptive benchmark. The
reported throughput is derived as `1000 / kg(ms)`.
