#ifndef FALCON_IMPL_OPT_CONFIG_H
#define FALCON_IMPL_OPT_CONFIG_H
#ifndef FALCON_IMPL_MIRROR_CACHE
#define FALCON_IMPL_MIRROR_CACHE 1
#endif
#ifndef FALCON_IMPL_VECTOR_NTT
#define FALCON_IMPL_VECTOR_NTT 1
#endif
#ifndef FALCON_IMPL_VECTOR_CDT
#if defined(__aarch64__)
#define FALCON_IMPL_VECTOR_CDT 1
#else
/* The original scalar scan is faster on the measured i7-13700. */
#define FALCON_IMPL_VECTOR_CDT 0
#endif
#endif
#ifndef FALCON_IMPL_SKIP_DEAD_OUTPUT
/* Exact, but not yet measured on the M5/i7 used for the selected defaults. */
#define FALCON_IMPL_SKIP_DEAD_OUTPUT 0
#endif
#if (FALCON_IMPL_MIRROR_CACHE != 0 && FALCON_IMPL_MIRROR_CACHE != 1) || \
    (FALCON_IMPL_VECTOR_NTT != 0 && FALCON_IMPL_VECTOR_NTT != 1) || \
    (FALCON_IMPL_VECTOR_CDT != 0 && FALCON_IMPL_VECTOR_CDT != 1) || \
    (FALCON_IMPL_SKIP_DEAD_OUTPUT != 0 && FALCON_IMPL_SKIP_DEAD_OUTPUT != 1)
#error "Falcon implementation switches must be 0 or 1"
#endif
#endif
