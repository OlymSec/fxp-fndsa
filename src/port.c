/* ANSA SDK platform layer; see antrag_port.h. */
#if defined(__linux__)
#define _GNU_SOURCE
#endif
#if defined(__APPLE__)
#define _DARWIN_C_SOURCE
#endif
#include "antrag_port.h"
#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#if defined(__linux__) || defined(__APPLE__)
#include <sys/random.h> /* getrandom (Linux), getentropy (macOS) */
#endif
#if defined(__APPLE__)
#include <dlfcn.h>
#include <sys/sysctl.h>
#endif

void antrag_port_cleanse(void *p, size_t len)
{
    if (!p || !len)
        return;
    volatile unsigned char *v = p;
#if (defined(__GNUC__) || defined(__clang__)) && UINTPTR_MAX > UINT32_MAX
    typedef uint64_t wipe_word __attribute__((__may_alias__));
    while (len && ((uintptr_t)v & 7U)) {
        *v++ = 0;
        len--;
    }
    volatile wipe_word *w = (volatile wipe_word *)(void *)v;
    while (len >= sizeof *w) {
        *w++ = 0;
        len -= sizeof *w;
    }
    v = (volatile unsigned char *)(void *)w;
#endif
    while (len--)
        *v++ = 0;
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : : "r"(p) : "memory");
#endif
}

void *antrag_port_malloc(size_t len)
{
    return malloc(len ? len : 1);
}

void *antrag_port_zalloc(size_t len)
{
    return calloc(1, len ? len : 1);
}

void antrag_port_free(void *p)
{
    free(p);
}

void antrag_port_clear_free(void *p, size_t len)
{
    if (!p)
        return;
    antrag_port_cleanse(p, len);
    free(p);
}

int antrag_port_run_once(pthread_once_t *once, void (*init)(void))
{
    return pthread_once(once, init) == 0;
}

typedef struct {
    void *buffer[ANTRAG_PORT_BUFFERS];
    size_t size[ANTRAG_PORT_BUFFERS];
} thread_buffers;

static pthread_key_t buffers_key;
static pthread_once_t buffers_once = PTHREAD_ONCE_INIT;
static int buffers_key_ok;
static _Thread_local thread_buffers *buffers_self;

/* Buffers of at least this size are mapped directly, so that releasing them
 * always returns their memory to the operating system. */
#define BUFFER_MAP_BYTES ((size_t)64 << 10)

static void *buffer_alloc(size_t size)
{
    if (size < BUFFER_MAP_BYTES)
        return calloc(1, size);
    void *p = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    return p == MAP_FAILED ? NULL : p; /* fresh anonymous pages are zero */
}

static size_t buffer_used(unsigned slot, void *buffer, size_t size);

/* Wipes the part of the buffer that may hold data, then frees it. */
static void buffer_free(unsigned slot, void *buffer, size_t size)
{
    if (!buffer)
        return;
    antrag_port_cleanse(buffer, buffer_used(slot, buffer, size));
    if (size < BUFFER_MAP_BYTES)
        free(buffer);
    else
        munmap(buffer, size);
}

static void buffers_destroy(void *pointer)
{
    thread_buffers *t = pointer;
    if (!t)
        return;
    buffers_self = NULL; /* runs on the exiting thread; later destructors must not see it */
    for (unsigned i = 0; i < ANTRAG_PORT_BUFFERS; i++)
        buffer_free(i, t->buffer[i], t->size[i]);
    antrag_port_clear_free(t, sizeof *t);
}

/* pthread key destructors do not run for the thread that calls exit() (or
 * returns from main), so that thread's buffers are wiped and freed here. */
static void buffers_at_exit(void)
{
    thread_buffers *t = buffers_self;
    if (!t)
        return;
    if (buffers_key_ok)
        pthread_setspecific(buffers_key, NULL);
    buffers_destroy(t);
}

/* The thread-exit destructor and the exit handler point into this library,
 * so it must stay loaded once they are registered. The Linux shared library
 * is linked with -z nodelete; on macOS the image pins itself. */
static void pin_image(void)
{
#if defined(__APPLE__)
    Dl_info info;
    if (dladdr((const void *)pin_image, &info) && info.dli_fname)
        dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE);
#endif
}

static void buffers_make_key(void)
{
    buffers_key_ok = pthread_key_create(&buffers_key, buffers_destroy) == 0;
    if (buffers_key_ok) {
        pin_image();
        atexit(buffers_at_exit);
    }
}

void *antrag_port_thread_buffer(unsigned slot, size_t size)
{
    if (slot >= ANTRAG_PORT_BUFFERS || !size)
        return NULL;
    thread_buffers *t = buffers_self;
    if (!t) {
        if (pthread_once(&buffers_once, buffers_make_key) != 0 || !buffers_key_ok)
            return NULL;
        t = calloc(1, sizeof *t);
        if (!t)
            return NULL;
        if (pthread_setspecific(buffers_key, t) != 0) {
            free(t);
            return NULL;
        }
        buffers_self = t;
    }
    if (!t->buffer[slot]) {
        t->buffer[slot] = buffer_alloc(size);
        if (!t->buffer[slot])
            return NULL;
        t->size[slot] = size;
    }
    return t->size[slot] == size ? t->buffer[slot] : NULL;
}

void *antrag_port_thread_buffer_peek(unsigned slot)
{
    thread_buffers *t = buffers_self;
    return t && slot < ANTRAG_PORT_BUFFERS ? t->buffer[slot] : NULL;
}

typedef struct {
    size_t top;        /* offset of the first free byte of data */
    size_t high_water; /* largest top so far: bytes that may hold data */
    unsigned char pad[48];
    unsigned char data[];
} scratch_region;

/* Only the scratch region's used prefix needs wiping; wiping the rest would
 * fault in pages that were never touched. */
static size_t buffer_used(unsigned slot, void *buffer, size_t size)
{
    if (slot != ANTRAG_PORT_BUFFER_SCRATCH)
        return size;
    scratch_region *r = buffer;
    return sizeof *r +
           (r->high_water < ANTRAG_PORT_SCRATCH_BYTES ? r->high_water : ANTRAG_PORT_SCRATCH_BYTES);
}

void *antrag_port_scratch_push(size_t size, int zero)
{
    scratch_region *r = antrag_port_thread_buffer(
        ANTRAG_PORT_BUFFER_SCRATCH, sizeof(scratch_region) + ANTRAG_PORT_SCRATCH_BYTES);
    if (!r)
        return NULL;
    size_t start = (r->top + 63u) & ~(size_t)63u;
    if (start > ANTRAG_PORT_SCRATCH_BYTES || size > ANTRAG_PORT_SCRATCH_BYTES - start)
        return NULL;
    unsigned char *block = r->data + start;
    if (zero)
        memset(block, 0, size);
    r->top = start + size;
    if (r->top > r->high_water)
        r->high_water = r->top;
    return block;
}

void antrag_port_scratch_pop(void *block, size_t size, int wipe)
{
    scratch_region *r = antrag_port_thread_buffer_peek(ANTRAG_PORT_BUFFER_SCRATCH);
    if (!r || !block)
        return;
    if (wipe)
        antrag_port_cleanse(block, size);
    r->top = (size_t)((unsigned char *)block - r->data); /* LIFO */
}

void antrag_port_thread_buffer_release(unsigned slot)
{
    thread_buffers *t = buffers_self;
    if (!t || slot >= ANTRAG_PORT_BUFFERS || !t->buffer[slot])
        return;
    buffer_free(slot, t->buffer[slot], t->size[slot]);
    t->buffer[slot] = NULL;
    t->size[slot] = 0;
}

#if !defined(__APPLE__) && !defined(__OpenBSD__) && !defined(__FreeBSD__) && !defined(__NetBSD__)
/* Fallback for systems without getentropy, and Linux kernels before getrandom. */
static int read_urandom(unsigned char *out, size_t len)
{
    unsigned char *start = out;
    int fd;
    do
        fd = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    while (fd < 0 && errno == EINTR);
    if (fd < 0)
        return 0;
    while (len) {
        ssize_t r = read(fd, out, len);
        if (r < 0 && errno == EINTR)
            continue;
        if (r <= 0) {
            close(fd);
            antrag_port_cleanse(start, (size_t)(out - start));
            return 0;
        }
        out += r;
        len -= (size_t)r;
    }
    close(fd);
    return 1;
}
#endif

int antrag_port_random(void *out, size_t len)
{
    unsigned char *p = out;
    if (!p && len)
        return 0;
    size_t total = len;
#if defined(__linux__)
    while (len) {
        ssize_t r = getrandom(p, len, 0);
        if (r < 0 && errno == EINTR)
            continue;
        if (r < 0 && errno == ENOSYS) {
            if (read_urandom(p, len))
                return 1;
            antrag_port_cleanse(out, total);
            return 0;
        }
        if (r <= 0) {
            antrag_port_cleanse(out, total);
            return 0;
        }
        p += r;
        len -= (size_t)r;
    }
    return 1;
#elif defined(__APPLE__) || defined(__OpenBSD__) || defined(__FreeBSD__) || defined(__NetBSD__)
    while (len) {
        size_t n = len < 256 ? len : 256; /* getentropy limit */
        if (getentropy(p, n) != 0) {
            antrag_port_cleanse(out, total);
            return 0;
        }
        p += n;
        len -= n;
    }
    return 1;
#else
    if (read_urandom(p, len))
        return 1;
    antrag_port_cleanse(out, total);
    return 0;
#endif
}

#if defined(__aarch64__) && (defined(__APPLE__) || defined(__linux__))
/* MSR DIT, #1, encoded so that no assembler extension is needed. */
#define DIT_SET ".inst 0xd503415f"
#if defined(__linux__)
#include <sys/auxv.h>
#ifndef HWCAP_DIT
#define HWCAP_DIT (1UL << 24)
#endif
static int dit_state; /* 0: not yet checked, 1: available, -1: not available */
static int dit_available(void)
{
    int state = __atomic_load_n(&dit_state, __ATOMIC_RELAXED);
    if (!state) {
        state = (getauxval(AT_HWCAP) & HWCAP_DIT) ? 1 : -1;
        __atomic_store_n(&dit_state, state, __ATOMIC_RELAXED);
    }
    return state > 0;
}
#else
#include <TargetConditionals.h>
#if TARGET_OS_OSX
/* Every Apple silicon Mac (M1 and later) implements DIT. */
static int dit_available(void)
{
    return 1;
}
#else
/* iOS, tvOS, watchOS: A11 and older processors lack DIT. */
static int dit_state; /* 0: not yet checked, 1: available, -1: not available */
static int dit_available(void)
{
    int state = __atomic_load_n(&dit_state, __ATOMIC_RELAXED);
    if (!state) {
        int enabled = 0;
        size_t size = sizeof enabled;
        if (sysctlbyname("hw.optional.arm.FEAT_DIT", &enabled, &size, NULL, 0))
            enabled = 0;
        state = enabled ? 1 : -1;
        __atomic_store_n(&dit_state, state, __ATOMIC_RELAXED);
    }
    return state > 0;
}
#endif
#endif

void antrag_port_dit_enable(void)
{
    if (dit_available())
        __asm__ volatile(DIT_SET ::: "memory");
}
#else
void antrag_port_dit_enable(void)
{
}
#endif
