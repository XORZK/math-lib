#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#pragma once
#include <cassert>
#include <cmath>
#include <iostream>

// a + bi
typedef struct complex {
	double a,b;

	complex();

	complex(const double A);

	complex(const double A, const double B);

	complex(const complex &c);

	~complex() {}

	double magnitude() const;

	complex conjugate() const;

	complex operator*=(const double scalar);

	complex operator/=(const double scalar);

	complex operator+=(const complex &c2);

	complex operator-=(const complex &c2);

	complex operator*=(const complex &c2);

	complex operator/=(const complex &c2);
} complex;

complex operator*(const double scalar, const complex &c1);

complex operator*(const complex &c1, const double scalar);

complex operator/(const double scalar, const complex &c1);

complex operator/(const complex &c1, const double scalar);

complex operator+(const complex &c1, const complex &c2);

complex operator-(const complex &c1, const complex &c2);

complex operator*(const complex &c1, const complex &c2);

complex operator/(const complex &c1, const complex &c2);

complex operator/(const double A, const complex &c);

bool operator==(const complex &c1, const complex &c2);

bool operator!=(const complex &c1, const complex &c2);

complex complex_sqrt(const complex &z);

complex complex_sqrt(const double A);

complex complex_pow(const complex &w, const complex &z);

complex complex_pow(const complex &w, const double b);

complex complex_pow(const double c, const complex &w);

complex complex_log(const complex &w);

std::ostream& operator<<(std::ostream &out, const complex &c);

#endif
