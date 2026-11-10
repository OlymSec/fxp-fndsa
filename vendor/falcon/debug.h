#include "fpr.h"

#ifdef DEBUG
#include <stdio.h>

extern char do_debug;

static inline void set_debug(char toggle) {
	do_debug = toggle;
}

static inline double approx(fpr x) {
#if FALCON_FXP
	fpr64 v = fpr_to_fpr64(x);
    double ret = ((double) *(int64_t*) &(v.v) / (double)((int64_t)  1 << 32));
    uint128_t s = x.v>>127;
    if (s & (ret >0) || !s & (ret < 0))
        printf("-APPROX ERROR (%016lx|%016lx\t%ld.)\t", (uint64_t)(x.v >> 64) , (uint64_t)x.v, (uint64_t)(x.v >> 64));
#elif FALCON_FPNATIVE
    double ret = x.v;
#elif FALCON_FPEMU
    double ret = *(double *)&x;
#endif
    return ret;
}

/* 
 * Debugging functions
*/
static inline void print_buf32(uint32_t *f, size_t len, char* label) {
	if(!do_debug)
		return;
	if (label) 
		printf("%s: ", label);
	for (size_t i = 0; i < len; i++)
		printf("%08x ", f[i]);
	printf("\n");
	fflush(stdout);
}

static inline void print_buf16(uint16_t *f, size_t len, char* label) {
	if(!do_debug)
		return;
	if (label)
		printf("%s: ", label);
	for (size_t i = 0; i < len; i++)
		printf("%04x ", (uint16_t)f[i]);
	printf("\n");
	fflush(stdout);
}

static inline void print_buf8(uint8_t *f, size_t len, char* label) {
	if(!do_debug)
		return;
	if (label)
		printf("%s: ", label);
	for (size_t i = 0; i < len; i++)
		printf("%02x ", (uint8_t)f[i]);
	printf("\n");
	fflush(stdout);
}

static inline void print_poly(int8_t *f, char *label) {
	if(!do_debug)
		return;
	print_buf8((uint8_t*)f, 512, label);
}

static inline void print_fp(fpr x) {
	if(!do_debug)
		return;
#if FALCON_FXP
    printf("%016lx|%016lx ", (uint64_t)(x.v >> 64) , (uint64_t)x.v);
#elif FALCON_FPNATIVE
    printf("%.15lf\n", x.v);
#elif FALCON_FPEMU
    printf("%.15lf\n", *(double *)&x);
#endif
}

static inline void print_fp_buf(fpr *fp, size_t len,  char* label) {
	if(!do_debug)
		return;
	if (label)
		printf("%s: ", label);
	for (size_t i = 0; i < len; i++) {
		print_fp(*(fp + i));
	}
	printf("\n");
	fflush(stdout);
}

static inline void print_fpapprox_buf(fpr *fp, size_t len,  char* label) {
	if(!do_debug)
		return;
	if (label)
		printf("%s: ", label);
	for (size_t i = 0; i < len; i++)
		printf("%.15lf\t", approx(*(fp + i)));
	printf("\n");
	fflush(stdout);
}

static inline void print_marker(char *s) {
	if(!do_debug)
		return;
	printf("#################################################\n");
	printf("#        %s             #\n", s);
	printf("#################################################\n");
}

#define print_debug(x) if (do_debug) {printf x;}

#else 

#define set_debug(x) 
#define approx(x)
#define print_buf32(x, y, z)
#define print_buf16(x, y, z)
#define print_buf8(x, y, z)
#define print_poly(x, y)
#define print_fp_buf(x, y, z)
#define print_fpapprox_buf(x, y, z)
#define print_marker(x)
#define print_debug(x)

#endif