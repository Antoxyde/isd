#ifndef _OP_339_
#define _OP_339_

#include <stdint.h>
#include <stddef.h>

void tests_mul339(void);
void tests_hash339(void);

uint64_t mul339(const uint64_t a, const uint64_t b);
uint64_t hash339(const uint32_t k, const void* buf, const size_t buflen);

#endif
