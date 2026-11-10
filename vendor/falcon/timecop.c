
#ifdef TIMECOP
#include <stddef.h>
#include "timecop/poison.h"
void
timecop_poison(void * addr, size_t len)
{
    poison(addr, len);
}

void
timecop_unpoison(void * addr, size_t len)
{
    unpoison(addr, len);
}

int
timecop_is_poisoned(void * addr, size_t len)
{
    return is_poisoned(addr, len);
}
#endif