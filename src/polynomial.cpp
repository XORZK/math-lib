#include "polynomial.hpp"

#define DEGREE(p) (p.N)
#define COEFF(p, N) (p.coeffs[N])

polynomial::polynomial() : N(0), coeffs(nullptr) {}

// degree N -> N+1 terms
polynomial::polynomial(size_t s) : N(s) {
	this->coeffs = static_cast<double*>(malloc(sizeof(double) * (N+1)));

	for (size_t k = 0; k <= N; k++) {
		this->coeffs[k] = 0.0;
	}
}

// Assumes c is array of size s+1
polynomial::polynomial(size_t s, double *c) : N(s) {
	this->coeffs = static_cast<double*>(malloc(sizeof(double) * (N+1)));

	for (size_t k = 0; k <= N; k++) {
		this->coeffs[k] = c[k];
	}
}

polynomial::polynomial(const polynomial &p) : N(p.N) {
	this->coeffs = static_cast<double*>(malloc(sizeof(double) * (N+1)));

	for (size_t k = 0; k <= N; k++) {
		this->coeffs[k] = p.coeffs[k];
	}
}

polynomial::~polynomial() {
	this->N = 0;
	//free(this->coeffs);
}

size_t polynomial::degree() const {
	return this->N;
}

void polynomial::coeff(const size_t deg, const double c) {
	assert(deg >= 0 && deg <= N);
	this->coeffs[deg] = c;
}

double polynomial::coeff(const size_t deg) const {
	assert(deg >= 0 && deg <= N);
	return this->coeffs[deg];
}

polynomial polynomial::derivative() const {
	polynomial d = polynomial(MAX(N-1, 0));

	// a_n * x^n becomes (a_n * n) * x^(n-1)
	for (size_t k = 1; k <= N; k++) {
		d.coeffs[k-1] = this->coeffs[k] * k;
	}

	return d;
}

double& polynomial::operator[](const size_t deg) {
	assert(deg >= 0 && deg <= N);
	return this->coeffs[deg];
}

// Evaluates P(x) where x is complex.
// TODO: Use Horner's Method
complex polynomial::operator()(const complex &x) {
	complex B = complex(this->coeffs[N], 0.0);

	for (int k = N-1; k >= 0; k--) {
		complex c = complex(this->coeffs[k], 0.0);
		B = c + B * x;
	}

	return B;
}

// Evaluates P(x) where x is real.
double polynomial::operator()(const double x) {
	double B = this->coeffs[N];

	for (int k = N-1; k >= 0; k--) {
		B = this->coeffs[k] + B * x; // Horner's method.
	}

	return B;
}

void polynomial::operator*=(const double scalar) {
	for (size_t k = 0; k <= N; k++)
		this->coeffs[k] *= scalar;
}

void polynomial::operator*=(polynomial &f) {
	polynomial h = ((*this) * f);
	if (h.N > this->N) {
		this->N = h.N;
		this->coeffs = static_cast<double*>(malloc(sizeof(double) * (N+1)));
	}

	for (size_t k = 0; k <= N; k++)
		this->coeffs[k] = h[k];
}

void polynomial::operator/=(const double scalar) {
	assert(scalar != 0);

	for (size_t k = 0; k <= N; k++)
		this->coeffs[k] /= scalar;
}

void polynomial::operator+=(polynomial &g) {
	double *tmp = static_cast<double*>(malloc(sizeof(double) * (N+1)));

	for (size_t k = 0; k <= MAX(this->N, g.N); k++) {
		double c = 0;
		if (k <= this->N)
			c += this->coeffs[k];
		if (k <= g.N)
			c += g[k];
		tmp[k] = c;
	}

	this->coeffs = static_cast<double*>(malloc(sizeof(double) * (N+1)));
	for (size_t k = 0; k <= N; k++)
		*(this->coeffs++) = *(tmp++);
}

void polynomial::operator-=(polynomial &g) {
	polynomial h = (g * -1);
	(*this) += h;
}

polynomial polynomial::trim(const double e) const {
	int M = this->N;

	while (ABS(this->coeffs[M]) < e && M > 0)
		M--;

	polynomial T(M);

	for (size_t k = 0; k <= M; k++)
		T[k] = this->coeffs[k];

	return T;
}

// (a_n)x^n + (a_(n-1))x^(n-1) + ... + a_0
// [a_n, a_(n-1),...,a_0]
std::ostream& operator<<(std::ostream& out, polynomial p) {
	int MNN = DEGREE(p);

	for (int k = DEGREE(p); k >= 0; k--) {
		if (p[k] != 0)
			MNN = k;
	}

	for (int k = DEGREE(p); k >= MNN; k--) {
		if (p[k] == 0.0)
			continue;

		out << std::format("{}({})(x^{})", (p[k] >= 0 ? "+" : "-"), ABS(p[k]), k);
	}

	return out;
}

polynomial power_func(size_t degree) {
	polynomial P(degree);
	P[degree] = 1;
	return P;
}

bool is_power_func(const polynomial &p) {
	polynomial T = p.trim();

	for (size_t k = 0; k < T.N; k++) {
		if (T[k] != 0)
			return false;
	}

	return true;
}

// returns the degree of the smallest degree, non-zero term
int min_degree(polynomial &p) {
	int M = DEGREE(p);

	for (int k = DEGREE(p); k >= 0; k--) {
		if (p[k] != 0)
			M = k;
	}

	return M;
}

bool operator==(const polynomial &f, const double value) {
	polynomial g = f.trim();
	return (g.N == 0 && g[0] == value);
}

polynomial operator*(const polynomial &f, const double scalar) {
	size_t N = DEGREE(f);

	polynomial g(N);

	for (size_t k = 0; k <= N; k++)
		g.coeff(k, f.coeff(k) * scalar);

	return g;
}

polynomial operator*(const double scalar, const polynomial &f) {
	return (f * scalar);
}

polynomial operator/(const polynomial &f, const double scalar) {
	assert(scalar != 0);

	return f * (1/scalar);
}

polynomial operator/(const double scalar, const polynomial &f) {
	assert(scalar != 0);

	return f * (1/scalar);
}

polynomial operator+(const polynomial &f, const polynomial &g) {
	size_t N = MAX(f.N, g.N);

	polynomial h(N);

	for (size_t k = 0; k <= N; k++) {
		double coeff = 0;

		if (k <= g.N)
			coeff += g.coeff(k);

		if (k <= f.N)
			coeff += f.coeff(k);

		h.coeff(k, coeff);
	}

	return h;
}

polynomial operator-(const polynomial &f, const polynomial &g) {
	return (f + (-1 * g));
}

// n^2 multiplications
// (n-1)^2 additions
polynomial schoolbook_algorithm(polynomial &A, polynomial &B) {
	polynomial C(A.N + B.N);

	for (size_t j = 0; j <= A.N; j++) {
		for (size_t k = 0; k <= B.N; k++)
			C[j+k] += (A[j] * B[k]);
	}

	return C;
}

// multiplies p by x^shift
polynomial shift_polynomial(polynomial p, const int shift) {
	int M = min_degree(p);
	// if dividing by x^|shift|, you want all non-zero coefficients
	// of p to have degree of at least |shift|
	if (shift < 0)
		assert(min_degree(p) >= ABS(shift));

	polynomial f(p.N + shift);

	for (int k = M; k <= p.N; k++)
		f[k+shift] = p[k];

	return f;
}

//					   upper       lower
// x^3 + x^2 + x + 1 = (x+1)*x^2 + (x+1)
// A = upper * x^k + lower
// basically, just take every term of degree at least k -> that is your upper.
// degree of lower is at most k-1.
// degree of upper is N-k.
void split_polynomial(polynomial &A, size_t k, polynomial *upper, polynomial *lower) {
	int N = DEGREE(A);

	*upper = polynomial(N-k);
	*lower = polynomial(k-1);

	for (size_t j = 0; j <= N; j++) {
		if (j < k)
			(*lower)[j] = A[j];
		else
			(*upper)[j-k] = A[j];
	}
}

// for polynomials of degree 2^n - 1
// for now, assume deg(A) = deg(B)
polynomial recursive_ka(polynomial &A, polynomial &B) {
	size_t N = MAX(A.N, B.N) + 1;

	if (N == 1)
		return schoolbook_algorithm(A, B);

	polynomial a_upper, a_lower, b_upper, b_lower;

	split_polynomial(A, N/2, &a_upper, &a_lower);
	split_polynomial(B, N/2, &b_upper, &b_lower);

	polynomial sum_a = a_upper + a_lower, sum_b = b_upper + b_lower;

	polynomial D0 = recursive_ka(a_lower, b_lower),
			   D1 = recursive_ka(a_upper, b_upper),
			   D01 = recursive_ka(sum_a, sum_b);

	return shift_polynomial(D1, N) + shift_polynomial(D01 - D0 - D1, N/2) + D0;
}

// deg(A) = deg(B) = 2
// A(x) = (a_2)*x^2 + (a_1)*x + a_0
// B(x) = (b_2)*x^2 + (b_1)*x + b_0
polynomial ka_quadratic(polynomial &A, polynomial &B) {
	polynomial C(4);

	double D0 = A[0] * B[0],
		   D1 = A[1] * B[1],
		   D2 = A[2] * B[2];

	double D01 = (A[0] + A[1])*(B[0] + B[1]),
		   D02 = (A[0] + A[2])*(B[0] + B[2]),
		   D12 = (A[1] + A[2])*(B[1] + B[2]);

	C[0] = D0;
	C[1] = (D01 - D0 - D1);
	C[2] = (D02 - D2 - D0 + D1);
	C[3] = (D12 - D1 - D2);
	C[4] = D2;

	return C;
}

// fast polynomial multiplication
// https://eprint.iacr.org/2006/224.pdf
// for now, assume deg(f) = deg(g) = N
polynomial ka_equal_degree(polynomial &f, polynomial &g) {
	assert(f.N == g.N);
	size_t N = f.N;

	if ((N+1) & N == 0)
		return recursive_ka(f, g);

	polynomial C(2*N);

	std::vector<double> a(N+1);
	std::vector<std::vector<double>> b;

	// index k stores D_k
	for (size_t k = 0; k <= N; k++) {
		std::vector<double> e(N+1);
		a[k] = (f[k] * g[k]);
		b.push_back(e);
	}

	// m from 1 -> 2(d+1)-3 = 2d-1
	// j + k = m
	// k > j >= 0
	for (int j = 0; j < b.size(); j++) {
		// (j, k) stores D_(j,k)
		for (int k = j+1; k < b[j].size(); k++) {
			b[j][k] = (f[j] + f[k])*(g[j] + g[k]);
		}
	}

	C[0] = a[0];
	C[2*N] = a[N];

	for (size_t i = 1; i < 2*N; i++) {
		double c = 0.0;
		for (int s = 0; s <= i/2+1 && i > s+s; s++) {
			int t = i-s;
			if (t >= N+1)
				continue;
			c += b[s][t] - (a[s] + a[t]);
		}

		if (i%2 == 0)
			c += a[i/2];

		C[i] = c;
	}

	return C;
}

polynomial operator*(polynomial &f, polynomial &g) {
	if (g.N > f.N)
		return (g*f);
	else if (f.N > g.N) {
		// basically just resizes the polynomial
		// and then shrinks it again
		polynomial T(f.N + g.N), h(f.N);
		for (size_t k = 0; k <= g.N; k++)
			h[k] = g.coeff(k);

		polynomial C = ka_equal_degree(f,h);

		for (size_t k = 0; k <= f.N + g.N; k++)
			T[k] = C[k];

		return T;
	}

	return ka_equal_degree(f, g);
}

// Divide P by monomial (x-x0)
// P(x) = Q(x)(x-x0) + R => Q(x) is returned
polynomial poly_divide(polynomial &P, const double x0) {
	size_t N = DEGREE(P);
	double *arr_b = static_cast<double*>(malloc(sizeof(double) * (N+1)));
	arr_b[N] = COEFF(P, N);

	for (int k = N-1; k >= 0; k--)
		arr_b[k] = COEFF(P, k) + arr_b[k+1] * x0;

	polynomial Q = polynomial(N-1);

	for (int k = N; k >= 1; k--)
		Q[k-1] = arr_b[k];

	free(arr_b);

	return Q;
}

// Divide polynomial P by monomial (x-x0) and load remainder into R
polynomial poly_divide(polynomial &P, const double x0, double *R) {
	*R = P(x0);
	return poly_divide(P, x0);
}

void poly_divide(polynomial &P, const double x0, polynomial *Q, double *R) {
	*Q = poly_divide(P, x0);
	*R = P(x0);
}

// Standard Euclidean Division algorithm
void poly_divide(polynomial &P, polynomial &f, polynomial *Q, polynomial *R) {
	polynomial divisor = f.trim(), curr = P.trim();
	*Q = polynomial((curr.N >= divisor.N ? curr.N - divisor.N : 0));

	size_t N = DEGREE(divisor);

	double leading = divisor[N];

	while (curr.N >= N && curr != 0) {
		double value = curr[curr.N]/leading;
		polynomial power = value * power_func(curr.N - N);
		(*Q)[curr.N - N] = value;
		curr = (curr - (power * divisor)).trim();
	}

	*R = polynomial(curr).trim();
}

// Solves ax^2 + bx + c = 0
void solve_quadratic(double a, double b, double c, complex *x0, complex *x1) {
	double D = b*b - 4*a*c;
	*x0 = (complex(-b, 0.0) + complex_sqrt(D))/(2*a);
	*x1 = (complex(-b, 0.0) - complex_sqrt(D))/(2*a);
}

// Uses Cardano's Method (https://proofwiki.org/wiki/Cardano's_Formula)
void solve_cubic(double a, double b, double c, double d,
				 complex *x0, complex *x1, complex *x2) {
	double Q = (3*a*c - b*b)/(9*a*a); // 2/9
	double R = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a); // -20/54 = -10/27
	double Q3 = Q*Q*Q;
	double R2 = R*R;
	double D = Q3 + R2;

	double S = R + std::sqrt(D);
	double T = R - std::sqrt(D);

	S = (S < 0 ? -std::pow(-S, 1.0/3.0) : std::pow(S, 1.0/3.0));
	T = (T < 0 ? -std::pow(-T, 1.0/3.0) : std::pow(T, 1.0/3.0));

	double X = S + T - (b/(3*a));
	double A = -(S+T)/2.0 - (b/(3*a));
	double B = std::sqrt(3.0/4.0) * (S-T);

	*x0 = complex(X, 0.0);
	*x1 = complex(A, B);
	*x2 = complex(A, -B);
}

void solve_cubic(polynomial &P, complex *x0, complex *x1, complex *x2) {
	assert(DEGREE(P) == 3);
	solve_cubic(COEFF(P, 3), COEFF(P, 2), COEFF(P, 1), COEFF(P, 0), x0, x1, x2);
}

// Recursively computes divided differences
double divided_differences_recurse(std::vector<double> &x, std::vector<double> &y, size_t k, size_t j) {
	assert(x.size() == y.size());

	if (k == j)
		return y[k];

	return (divided_differences_recurse(x, y, k+1, j) - divided_differences_recurse(x, y, k, j-1))/(x[k+j] - x[k]);
}

double divided_differences_recurse(std::vector<double> &x, std::vector<double> &y) {
	assert(x.size() == y.size());

	return divided_differences_recurse(x, y, 0, x.size()-1);
}

// Returns 2-D vector diffs
// diffs[j][k] = [y_k, ..., y_(k+j)]
std::vector<std::vector<double>> divided_differences(std::vector<double> &x, std::vector<double> &y) {
	std::vector<double>::size_type N = x.size();
	std::vector<std::vector<double>> diffs;

	for (size_t k = 0; k < N; k++) {
		std::vector<double> tmp(N-k);
		diffs.push_back(tmp);
	}

	// N(N+1)/2 => O(N^2)
	for (size_t j = 0; j < N; j++) {
		// k+j < N => j < N-k
		for (size_t k = 0; k < N-j; k++) {
			assert(j == 0 || x[k+j] != x[k]);
			if (j == 0)
				diffs[0][k] = y[k];
			else
				diffs[j][k] = (diffs[j-1][k+1] - diffs[j-1][k])/(x[k+j] - x[k]); // indexes at (j-1, k) and (j-1, k+1) exist
		}
	}

	return diffs;
}

std::vector<std::vector<double>> divided_differences(std::vector<std::pair<double, double>> &points) {
	std::vector<double>::size_type N = points.size();
	std::vector<std::vector<double>> diffs;

	for (size_t k = 0; k < N; k++) {
		std::vector<double> tmp(N-k);
		diffs.push_back(tmp);
	}

	for (size_t j = 0; j < N; j++) {
		for (size_t k = 0; k < N-j; k++) {
			assert(j == 0 || points[k+j].first != points[k].first);
			if (j == 0)
				diffs[0][k] = points[k].second;
			else
				diffs[j][k] = (diffs[j-1][k+1] - diffs[j-1][k])/(points[k+j].first - points[k].first);
		}
	}

	return diffs;
}

// two values a0, b0 such that f(a0)f(b0) < 0 (opposite signs)
// if f is continuous, a solution is guaranteed on (a0, b0)
// stop when |a-b| < e.
double dekkers_method(polynomial &f, const double a0, const double b0, const double e) {
	assert((f(a0)*f(b0)) < 0);

	// a: a_k, b: b_k, c: b_(k-1)
	double a = a0, b = b0, c = a0;

	while (ABS(a-b) >= e) {
		if (ABS(f(a)) < ABS(f(b)))
			std::swap(a, b);

		double next_a, next_b, next_yb;
		double s, m = (a+b)/2;
		double ya = f(a), yb = f(b), yc = f(c);

		if (ya != yb)
			s = b-yb*(b-c)/(yb-yc);
		else
			s = (a+b)/2;

		if ((m > b && b < s && s < m) || (m < b && m < s && s < b))
			next_b = s;
		else
			next_b = m;

		next_yb = f(next_b);
		if ((ya*next_yb) < 0)
			next_a = a;
		else
			next_a = b;

		c = b;
		a = next_a;
		b = next_b;
	}

	return b;
}

// assume f(a0) < 0 and f(b0) > 0
double bisection_method(polynomial &f, const double a0, const double b0, const double e) {
	assert((f(a0) * f(b0)) < 0);

	double a = a0, b = b0;

	if (f(a0) > 0)
		std::swap(a, b);

	while (ABS(a-b) >= e) {
		double m = (a+b)/2.0;
		double ym = f(m);

		if (ym > 0)
			b = m;
		else
			a = m;

		// does |f(m)| < e make more sense.
		if (ABS(ym) < e)
			return m;
	}

	return b;
}

// stop when |a-b| < e
// https://en.wikipedia.org/wiki/Brent%27s_method
double brents_method(polynomial &f, const double a0, const double b0, const double tol, const double e) {
	assert((f(a0) * f(b0)) < 0);

	bool prev_bisect = false;

	// a: a_k, b: b_k, c: b_(k-1), d: b_(k-2)
	double a = a0, b = b0, c = a0, d, s;

	while (ABS(b-a) >= e) {
		double ya = f(a), yb = f(b), yc = f(c), ys;
		double bounds = (3*a + b)/4;

		if (ya != yb && yb != yc && ya != yc)
			// IQI: https://en.wikipedia.org/wiki/Inverse_quadratic_interpolation
			s = a*(yb*yc)/((ya-yb)*(ya-yc)) + b*(ya*yc)/((yb-ya)*(yb-yc)) + c*(ya*yb)/((yc-ya)*(yc-yb));
		else
			// Secant Method: https://en.wikipedia.org/wiki/Secant_method
			s = b - yb*(b-a)/(yb-ya);

		if (prev_bisect && ((ABS(tol) >= ABS(b-c)) || (ABS(s-b) >= (0.5*ABS(b-c)))) ||
			!prev_bisect && ((ABS(tol) >= ABS(c-d)) || (ABS(s-b) >= (0.5*ABS(c-d)))) ||
			((bounds < b && !(bounds < s && s < b)) || (b < bounds && !(b < s && s < bounds)))) {
			s = (a+b)/2;
			prev_bisect = true;
		} else
			prev_bisect = false;

		ys = f(s);

		d = c;
		c = b;

		if ((ya*ys) < 0)
			b = s;
		else
			a = s;

		if (ABS(f(a)) < ABS(f(b)))
			std::swap(a, b);

		if (ys == 0)
			break;
	}

	return b;
}

// Attempts at generating a Sturm chain
// https://web.math.ucsb.edu/~padraic/mathcamp_2013/root_find_alg/Mathcamp_2013_Root-Finding_Algorithms_Day_2.pdf
std::vector<polynomial> sturm_chain(polynomial &p) {
	std::vector<polynomial> chain;

	chain.push_back(p);
	chain.push_back(p.derivative());

	while (true) {
		polynomial Q, R;
		polynomial p0 = chain[chain.size()-1],
				   p1 = chain[chain.size()-2];
		poly_divide(p1, p0, &Q, &R);

		chain.push_back(-1 * R);

		if (DEGREE(R) == 0)
			break;
	}

	return chain;
}

size_t sign_changes(std::vector<polynomial> &chain, const double value) {
	size_t count = 0;
	std::vector<double> values;
	for (polynomial &p : chain) {
		double c = p(value);
		if (c != 0)
			values.push_back(c);
	}

	// counts sign variations
	// how many changes are there from + to - or - to +
	for (size_t k = 0; k < values.size() - 1; k++) {
		if ((values[k]*values[k+1]) < 0)
			count++;
	}

	return count;
}

size_t sign_changes(polynomial &p, const double value) {
	std::vector<polynomial> chain = sturm_chain(p);
	return sign_changes(chain, value);
}

// More efficient (Sturm chain only needs to be computed once).
size_t distinct_roots(std::vector<polynomial> &chain, const double a, const double b) {
	if (a > b)
		return distinct_roots(chain, b, a);

	size_t Va = sign_changes(chain, a),
		   Vb = sign_changes(chain, b);

	return (Va >= Vb ? Va - Vb : 0);
}

// Uses Sturm's theorem to count the number of distinct real roots of p
// in half-open interval [a, b) (a < b)
size_t distinct_roots(polynomial &p, const double a, const double b) {
	if (a > b)
		return distinct_roots(p, b, a);

	std::vector<polynomial> chain = sturm_chain(p);

	return distinct_roots(chain, a, b);
}

size_t distinct_roots(std::vector<polynomial> &chain) {
	std::vector<double> neg, pos;

	for (polynomial &f : chain) {
		double leading = f[f.N];
		// negative inf
		double n_inf = (leading * (f.N % 2 == 0) ? 1 : -1);

		if (!f.N)
			n_inf = leading;

		// positive inf
		double p_inf = (leading);
		neg.push_back(n_inf > 0 ? 1 : -1);
		pos.push_back(p_inf > 0 ? 1 : -1);
	}

	// compute V_a, V_b
	// neg.size() = pos.size()
	size_t V_a = 0, V_b = 0;
	for (size_t k = 0; k < neg.size()-1; k++) {
		if ((neg[k] * neg[k+1]) < 0)
			V_a++;
		if ((pos[k] * pos[k+1]) < 0)
			V_b++;
	}

	return (V_a >= V_b ? V_a - V_b : 0);
}

// computes the total amount of distinct roots of p
// on (-inf, inf)
size_t distinct_roots(polynomial &p) {
	std::vector<polynomial> chain = sturm_chain(p);
	return distinct_roots(chain);
}

void find_roots(polynomial &p, std::vector<polynomial> *chain, std::vector<double> *roots, const double a, const double b) {
	if (b <= a) {
		if (b == a && p(a) == 0)
			(*roots).push_back(a);
		return;
	}

	size_t count = distinct_roots(*chain, a, b);

	if (count == 0)
		return;

	if (count == 1 && ((p(b) * p(a)) < 0)) {
		(*roots).push_back(brents_method(p, a, b));
		return;
	}

	double m = (a+b)/2;

	find_roots(p, chain, roots, a, m);
	find_roots(p, chain, roots, m, b);
}

// Searches for roots a combination of Brent's method & Sturm's theorem.
//
// The idea is to perform some sort of divide and conquer algorithm.
// If we are searching in interval [a,b) (a < b), at each step compute
// some m such that m = (a+b)/2 and then search the interval [a,m] and
// [m, b] respectively. We can use Sturm's theorem to narrow down on the
// intervals containing roots. For example, if [a,b) contains n distinct
// roots, then [a,m) and [m,b) both contain at most n distinct roots.
// Also, if m is not a root, then count([a,m)) + count([m, b)) = n.
// If count([a,m)) = 0 or count([m,b)) = 0, we can simply discard
// that range for searching. Also, if p(a) * p(m) < 0 or p(m) * p(b) < 0,
// then we have two values that can be used for root searching by
// Brent's method. If not, we can just keep closing down the interval
// until [a,m) is so small or [m, b) is so small that a root must be at
// x = m.
std::vector<double> find_roots(polynomial &p, const double min_x, const double max_x, size_t max_its) {
	size_t its = 0;
	std::vector<double> roots(0);
	std::vector<polynomial> chain = sturm_chain(p);

	double total_roots = distinct_roots(chain);
	double a = min_x, b = max_x;

	if (a > b)
		std::swap(a,b);

	while (distinct_roots(chain,a,b) != total_roots && its <= max_its) {
		a -= 10000;
		b += 10000;
		its++;
	}

	find_roots(p, &chain, &roots, a, b);
	return roots;
}

// this assumes p is unimodal on [a0, b0]
double gss(polynomial &p, const double a0, const double b0, const bool min, const double tol) {
	double r = (-1 + std::sqrt(5))/2;

	double a = a0, b = b0;
	double x1, x2;
	double y1, y2;

	x1 = a + (1-r)*(b-a);
	x2 = a + (r)*(b-a);

	y1 = p(x1);
	y2 = p(x2);

	while (ABS(b-a) > tol) {
		if ((min && y1 > y2) || (!min && y1 < y2)) {
			a = x1;
			x1 = x2;
			y1 = y2;
			x2 = a + (r)*(b-a);
			y2 = p(x2);
		} else {
			b = x2;
			x2 = x1;
			y2 = y1;
			x1 = a + (1-r)*(b-a);
			y1 = p(x1);
		}
	}

	return (a+b)/2;
}
