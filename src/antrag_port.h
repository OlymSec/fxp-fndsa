/* ANSA SDK platform layer: the only operating-system services the Antrag
 * library needs (memory wiping, allocation, one-time initialization and
 * secure randomness). The research sources used OpenSSL for these; the
 * importer rewrites them to call this layer, so the SDK links against no
 * third-party library and uses no OpenSSL header or name. All functions are
 * hidden and localized in the built libraries. */
#ifndef ANTRAG_PORT_H
#define ANTRAG_PORT_H

#include <pthread.h>
#include <stddef.h>

/* Overwrites len bytes with zeros; the compiler cannot elide the stores. */
void antrag_port_cleanse(void *p, size_t len);
void *antrag_port_malloc(size_t len);
void *antrag_port_zalloc(size_t len);
void antrag_port_free(void *p);
/* Wipes len bytes, then frees. Accepts NULL. */
void antrag_port_clear_free(void *p, size_t len);
/* Runs init exactly once per flag. Returns 1 on success, 0 on failure. */
int antrag_port_run_once(pthread_once_t *once, void (*init)(void));
/* Fills out with len bytes from the operating system's CSPRNG. Returns 1 on
 * success and 0 on failure; it never returns partially filled output as
 * success and clears the output buffer on failure. */
int antrag_port_random(void *out, size_t len);

/* Data-independent timing for operations on secret keys. On AArch64,
 * antrag_port_dit_enable sets PSTATE.DIT, under which the architecture
 * guarantees data-independent timing for data-processing instructions (and
 * Apple M3 and later disable the data-memory-dependent prefetcher). Every
 * secret-dependent API call sets it; DIT is per-thread state and is left
 * set. macOS sets it unconditionally, since every Apple silicon Mac
 * implements DIT; elsewhere CPU support is checked once per process and
 * cached, because the instruction is undefined before ARMv8.4 (Cortex-A72,
 * Neoverse N1, A11, ...). On other targets it does nothing. */
void antrag_port_dit_enable(void);

/* Per-thread working buffers. Large signing and key-generation workspaces
 * used to be thread-local variables, which glibc allocates and zero-fills in
 * every thread of the process. Instead, a thread gets a buffer only when it
 * first needs it: the buffer is zero-initialized, persists for the thread's
 * lifetime (exactly like the thread-local variable it replaces), and is wiped
 * and freed when the thread exits (for the thread that calls exit(), by an
 * atexit handler). Each slot always has the same size. */
enum { ANTRAG_PORT_BUFFER_SCRATCH, ANTRAG_PORT_BUFFERS };
/* Returns this thread's buffer for slot, allocating it on first use, or NULL
 * if allocation fails. */
void *antrag_port_thread_buffer(unsigned slot, size_t size);
/* Returns this thread's buffer for slot if it exists, without allocating. */
void *antrag_port_thread_buffer_peek(unsigned slot);
/* Wipes and frees this thread's buffer for slot now (it is allocated again
 * on next use). The buffer must not be in use. */
void antrag_port_thread_buffer_release(unsigned slot);

/* Per-thread scratch region for arrays that used to live on the call stack
 * (patches/0002-stack.patch). It is ordinary heap memory, reserved once per
 * thread on first use and reused afterwards, so a call costs no allocator
 * work; untouched pages are not resident. Blocks are 64-byte aligned and
 * strictly LIFO: every push is matched by a pop of the same block in reverse
 * order. zero clears the block (for zero-initialized arrays); wipe
 * overwrites it on pop (for secret data). Push returns NULL if the region is
 * exhausted. The region is wiped and freed at thread exit. */
#define ANTRAG_PORT_SCRATCH_BYTES ((size_t)1 << 20)
void *antrag_port_scratch_push(size_t size, int zero);
void antrag_port_scratch_pop(void *block, size_t size, int wipe);

#endif
