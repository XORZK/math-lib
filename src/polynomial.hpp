#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#pragma once
#include <iostream>
#include <format>
#include <string>
#include <string_view>
#include <vector>
#include "complex.hpp"
#include "MACROS.hpp"

typedef struct polynomial {
	size_t N;
	double *coeffs;

	polynomial();

	polynomial(size_t s);

	polynomial(size_t s, double *c);

	polynomial(const polynomial &p);

	~polynomial();

	size_t degree() const;

	void coeff(const size_t deg, const double c);

	double coeff(const size_t deg) const;

	polynomial derivative() const;

	double& operator[](const size_t deg);

	complex operator()(const complex &x);

	double operator()(const double x);

	// removes any leading coefficients if
	// |ak| < e
	// 0x^3 + 2x^2 + x + 1 -> 2x^2 + x + 1
	polynomial trim(const double e=10e-10) const;
} polynomial;

std::ostream& operator<<(std::ostream& out, polynomial p);

polynomial power_func(size_t degree);

int min_degree(polynomial &p);

bool operator==(const polynomial &f, const double value);

polynomial operator*(const polynomial &f, const double scalar);

polynomial operator*(const double scalar, const polynomial &f);

polynomial operator/(const polynomial &f, const double scalar);

polynomial operator/(const double scalar, const polynomial &f);

polynomial operator+(const polynomial &f, const polynomial &g);

polynomial operator-(const polynomial &f, const polynomial &g);

polynomial schoolbook_algorithm(polynomial &A, polynomial &B);

void split_polynomial(polynomial &A, size_t k, polynomial *upper, polynomial *lower);

polynomial recursive_ka(polynomial &A, polynomial &B);

polynomial ka_quadratic(polynomial &A, polynomial &B);

polynomial ka_equal_degree(polynomial &f, polynomial &g);

polynomial operator*(polynomial &f, polynomial &g);

// Divide polynomial P by monomial (x-x0)
polynomial poly_divide(polynomial &P, const double x0);

// Divide polynomial P by monomial (x-x0) and load remainder into R
polynomial poly_divide(polynomial &P, const double x0, double *R);

void poly_divide(polynomial &P, const double x0, polynomial *Q, double *R);

// Divide polynomial P by polynomial f
void poly_divide(polynomial &P, polynomial &f, polynomial *Q, polynomial *R);

// Solves ax^2 + bx + c = 0
void solve_quadratic(double a, double b, double c, complex *x0, complex *x1);

void solve_cubic(double a, double b, double c, double d, complex *x0, complex *x1, complex *x2);

void solve_cubic(polynomial &P, complex *x0, complex *x1, complex *x2);

double divided_differences_recurse(std::vector<double> &x, std::vector<double> &y, size_t k, size_t j);

double divided_differences_recurse(std::vector<double> &x, std::vector<double> &y);

std::vector<std::vector<double>> divided_differences(std::vector<double> &x, std::vector<double> &y);

std::vector<std::vector<double>> divided_differences(std::vector<std::pair<double, double>> &points);

double bisection_method(polynomial &f, const double a0, const double b0, const double e=10e-7);

double dekkers_method(polynomial &f, const double a0, const double b0, const double e=10e-7);

double brents_method(polynomial &f, const double a0, const double b0, const double tol = 10e-4, const double e=10e-7);

std::vector<polynomial> sturm_chain(polynomial &p);

size_t sign_changes(std::vector<polynomial> &chain, const double value);

size_t sign_changes(polynomial &p, const double value);

size_t distinct_roots(std::vector<polynomial> &chain, const double a, const double b);

size_t distinct_roots(polynomial &p, const double a, const double b);

size_t distinct_roots(std::vector<polynomial> &chain);

size_t distinct_roots(polynomial &p);

// searches for roots
std::vector<double> find_roots(polynomial &p, const double min_x, const double max_x, size_t max_its = 10);

#endif
