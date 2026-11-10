#ifndef FALCON_RECONSTRUCTION_H
#define FALCON_RECONSTRUCTION_H
#include <stddef.h>
size_t study_prepared_size(void);
/* The expanded key must remain immutable while this cache is used. */
int study_prepare_expanded(void *cache,const void *expanded_key);
void study_select_prepared(const void *cache);
void study_clear_temporary(void (*clear)(void *,size_t));
#endif
