#ifndef RATIONAL_HPP
#define RATIONAL_HPP

#pragma once
#include <iostream>
#include <stdint.h>

typedef struct rational {
	int64_t P;
	uint64_t Q;

	rational();

	~rational();
} rational;

#endif
