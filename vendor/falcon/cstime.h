#ifndef FALCON_CSTIME_H__
#define FALCON_CSTIME_H__

/* 
 * Implementation of small gadgets used for constant time operations.
*/

#include <stdint.h>

/* 
 * Constant time conditional move: if b is true, we move v into *r
*/
static inline void cmov_uint64(uint64_t *r, uint64_t v, uint64_t b)
{
  b = -b;
  *r ^= b & ((*r) ^ v);
}


static inline void cmov_int(int64_t *r, uint64_t v, uint64_t b)
{
  b = -b;
  *r ^= b & ((*r) ^ v);
}



#endif