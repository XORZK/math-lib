#include <iostream>
#include "complex.hpp"
#include "polynomial.hpp"
#include "newton_poly.hpp"

int main(int argc, char* argv[]) {
	/*
	complex x0, x1, x2;

	polynomial p = polynomial(3);

	p[0] = 4;
	p[1] = 3;
	p[2] = 2;
	p[3] = 1;

	//solve_cubic(1.0, 2.0, 3.0, 4.0, &x0, &x1, &x2);
	solve_cubic(p, &x0, &x1, &x2);

	std::cout << x0 << "\n";
	std::cout << x1 << "\n";
	std::cout << x2 << "\n";

	std::vector<double> x = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
	std::vector<double> y = { 1.0, 18.0, 52.0, 102.0, 201.0, 400.0 };

	std::vector<std::vector<double>> a = divided_differences(x, y);

	for (std::vector<double> b : a) {
		for (double c : b) {
			std::cout << c << " ";
		}
		std::cout << "\n";
	}

	newton_poly p;

	p.add_point(1.0, 1.0);
	p.add_point(2.0, 18.0);
	p.add_point(3.0, 52.0);
	p.add_point(4.0, 102.0);
	p.add_point(5.0, 201.0);
	p.add_point(6.0, 400.0);

	std::cout << p << "\n";
	std::cout << p(7) << "\n";

	polynomial A(3), B(3);
	A[0] = 1;
	A[1] = 2;
	A[2] = 5;
	A[3] = 3;

	B[0] = 2;
	B[1] = 4;
	B[2] = 6;
	B[3] = 8;

	std::cout << A << "\n";
	std::cout << bisection_method(A, 0, -2) << "\n";
	std::cout << dekkers_method(A, 0, -2) << "\n";
	std::cout << brents_method(A, 0, -2) << "\n";

	polynomial u1, l1;
	polynomial u2, l2;

	split_polynomial(A, 2, &u1, &l1);
	split_polynomial(u1, 1, &u2, &l2);
	std::cout << recursive_ka(A, B) << "\n";

	polynomial q1(2), q2(2);
	q1[0] = 1;
	q1[1] = 4;
	q1[2] = 5;

	q2[0] = 2;
	q2[1] = 6;
	q2[2] = 10;

	std::cout << ka_equal_degree(A, B) << "\n";*/
	polynomial p(4);
	p[0] = -1;
	p[1] = -1;
	p[3] = 1;
	p[4] = 1;

	std::vector<polynomial> chain = sturm_chain(p);

	for (polynomial &p : chain)
		std::cout << p << "\n";

	std::cout << distinct_roots(p, -10000, 10000) << "\n";
	std::cout << distinct_roots(p) << "\n";
	std::vector<double> roots = find_roots(p, -100, 100);

	for (double &b : roots) {
		std::cout << std::format("Root at: x = {}", b) << "\n";
	}
}
