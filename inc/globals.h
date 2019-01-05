#pragma once

#include <stdint.h>

#define CAPACITY    4096
#define NHASHES     12
#define HASHSIZE    512     // HASHSIZE % 32 must be 0


typedef uint16_t count_dtype;
#define MAXCOUNT    UINT16_MAX
