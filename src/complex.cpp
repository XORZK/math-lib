#include "complex.hpp"
#include "MACROS.hpp"

#define REAL(c) (c.a)
#define IMAG(c) (c.b)
#define COMPLEX_ARG(c) (c.a == 0.0 && c.b == 0.0 ? 0 : atan2(c.b, c.a))

complex::complex() : a(0), b(0) {}

complex::complex(const double A) : a(A), b(A) {}

complex::complex(const double A, const double B) : a(A), b(B) {}

complex::complex(const complex &c) : a(c.a), b(c.b) {}

double complex::magnitude() const {
	return std::sqrt(this->a*this->a + this->b*this->b);
}

// a+bi --> a-bi
complex complex::conjugate() const {
	return complex(this->a, -this->b);
}

complex complex::operator*=(const double scalar) {
	this->a *= scalar;
	this->b *= scalar;

	return *this;
}

complex complex::operator/=(const double scalar) {
	assert(scalar != 0);
	this->a /= scalar;
	this->b /= scalar;

	return *this;
}

// (a + bi) + (c + di) = (a+c) + (b+d)i
complex complex::operator+=(const complex &c2) {
	this->a += c2.a;
	this->b += c2.b;

	return *this;
}

// (a + bi) - (c + di) = (a-c) + (b-d)i
complex complex::operator-=(const complex &c2) {
	this->a -= c2.a;
	this->b -= c2.b;

	return *this;
}

// (a+bi) * (c+di) = (ac - bd) + (ad + bc)i
complex complex::operator*=(const complex &c2) {
	double t_a = (this->a * c2.a - this->b * c2.b),
		   t_b = (this->a * c2.b + this->b * c2.a);
	this->a = t_a;
	this->b = t_b;

	return *this;
}

// w/z
complex complex::operator/=(const complex &c2) {
	complex c1 = complex(*this);
	complex tmp = c1 * c2;
	this->a = tmp.a;
	this->b = tmp.b;

	return *this;
}

// N * (a+bi) = (aN) + (bN)i
complex operator*(const double scalar, const complex &c1) {
	return complex(c1.a * scalar, c1.b * scalar);
}

// (a+bi) * N = (aN) + (bN)i
complex operator*(const complex &c1, const double scalar) {
	return complex(c1.a * scalar, c1.b * scalar);
}

// (a+bi)/N = (a/N) + (b/N)i
complex operator/(const complex &c1, const double scalar) {
	assert(scalar != 0);
	return complex(c1.a / scalar, c1.b / scalar);
}

// (a + bi) + (c + di) = (a+c) + (b+d)i
complex operator+(const complex &c1, const complex &c2) {
	return complex(c1.a + c2.a, c1.b + c2.b);
}

// (a + bi) - (c + di) = (a-c) + (b-d)i
complex operator-(const complex &c1, const complex &c2) {
	return complex(c1.a - c2.a, c1.b - c2.b);
}

// (c1.a + c1.bi) * (c2.a + c2.bi)
// = (c1.a * c2.a - c1.b * c2.b) + (c1.a * c2.b + c1.b + c2.a)i
complex operator*(const complex &c1, const complex &c2) {
	return complex(c1.a * c2.a - c1.b * c2.b,
				   c1.a * c2.b + c1.b * c2.a);
}

// w/z
complex operator/(const complex &c1, const complex &c2) {
	return (c1 * c2.conjugate())/(std::pow(c2.magnitude(), 2));
}

// 1/z
complex operator/(const double A, const complex &c) {
	return A * (c.conjugate()/(std::pow(c.magnitude(), 2)));
}

bool operator==(const complex &c1, const complex &c2) {
	return (c1.a == c2.a && c1.b == c2.b);
}

bool operator!=(const complex &c1, const complex &c2) {
	return !(c1.a == c2.a && c1.b == c2.b);
}

complex complex_sqrt(const complex &z) {
	assert(z.b != 0);

	double R = z.magnitude();
	double A = std::sqrt((z.a + R)/2),
		   B = SIGN(z.b) * std::sqrt((-z.a + R)/2);
	return complex(A, B);
}

complex complex_sqrt(const double A) {
	if (A >= 0)
		return complex(std::sqrt(A), 0);
	return complex(0, std::sqrt(-A));
}

// w^z
complex complex_pow(const complex &w, const complex &z) {
	double theta = COMPLEX_ARG(w);
	double ln_r = std::log(w.magnitude());
	double M = std::exp(z.a * ln_r - z.b * theta);
	double beta = (z.b * ln_r + theta * z.a);

	return complex(M * cos(beta), M * sin(beta));
}

// w^b
complex complex_pow(const complex &w, const double b) {
	double theta = COMPLEX_ARG(w);
	double ln_r = std::log(w.magnitude());
	double M = std::exp(ln_r * b);
	double beta = theta * b;

	return complex(M * cos(beta), M * sin(beta));
}

// c^w
complex complex_pow(const double c, const complex &w) {
	double ln_c = std::log(c);
	double theta = ln_c * w.b;

	return complex(std::exp(ln_c * w.a) * cos(theta),
				   std::exp(ln_c * w.a) * sin(theta));
}

// ln(z)
complex complex_log(const complex &w) {
	return complex(std::log(w.magnitude()), COMPLEX_ARG(w));
}

// ln(z) / ln(base)
complex complex_log(const complex &w, const double base) {
	return complex_log(w)/std::log(base);
}

std::ostream& operator<<(std::ostream& out, const complex& c) {
	return (out << std::to_string(c.a) + (c.b > 0 ? "+" : "-") + std::to_string(ABS(c.b)) + "i");
}
