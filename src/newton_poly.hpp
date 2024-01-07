#ifndef NEWTON_POLY_HPP
#define NEWTON_POLY_HPP

#pragma once
#include <algorithm>
#include <utility>
#include <vector>
#include "complex.hpp"
#include "MACROS.hpp"
#include "polynomial.hpp"

// A polynomial which is represented through Newton's Divided Difference representation.
typedef struct newton_poly {
	// N is the number of points, not the degree of the polynomial.
	size_t N;
	// 2d array is not the most efficient way to store diffs.
	// 1d array can be used (althought it is more complicated).
	// Using an 1d array would be significantly less
	// efficient when adding new points.
	std::vector<std::vector<double>> diffs;
	std::vector<std::pair<double, double>> points;

	newton_poly();

	newton_poly(const newton_poly &p);

	~newton_poly();

	complex basis(const complex &x, const size_t j);

	double basis(const double x, const size_t j);

	size_t degree() const;

	void add_point(const double x0, const double y0);

	void remove_point(const double x0, const double y0);

	complex operator()(const complex &x);

	double operator()(const double x);
} newton_poly;

std::ostream& operator<<(std::ostream& out, newton_poly &p);

#endif
