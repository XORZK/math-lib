#include "newton_poly.hpp"

newton_poly::newton_poly() : N(0) {}

newton_poly::newton_poly(const newton_poly &p) : N(p.N) {
	this->points = p.points;
	this->diffs = p.diffs;
}

newton_poly::~newton_poly() {}

complex newton_poly::basis(const complex &x, const size_t j) {
	complex val = complex(1.0, 0.0);

	for (size_t k = 0; k < j; k++) {
		val *= (x - complex(points[k].first, 0.0));
	}

	return val;
}

double newton_poly::basis(const double x, const size_t j) {
	double val = 1.0f;

	for (size_t k = 0; k < j; k++)
		val *= (x - points[k].first);

	return val;
}

size_t newton_poly::degree() const {
	return MAX(0, this->points.size() - 1);
}

void newton_poly::add_point(const double x0, const double y0) {
	// TODO: want to make sure that x0 is not already in points.
	//assert();

	this->N++;
	this->points.push_back(std::make_pair(x0,y0));

	std::vector<double> last_rank(0);
	this->diffs.push_back(last_rank);

	// recompute divided differences.
	// [y0] [y1] [y2]		[y0] [y1] [y2] *[y3]
	// [y0,y1][y1,y2] =>    [y0,y1] [y1,y2] *[y2,y3]
	// [y1,y2,y3]	        [y1,y2,y3] *[y2,y3,y4]
	//					   *[y1,y2,y3,y4]
	// only one new difference needs to be computed per row.
	for (size_t j = 0; j < N; j++) {
		size_t k = N - 1 - j;
		if (j == 0)
			diffs[j].push_back(y0);
		else
			diffs[j].push_back((diffs[j-1][k+1] - diffs[j-1][k])/(points[k+j].first-points[k].first));
	}
}

void newton_poly::remove_point(const double x0, const double y0) {
	std::pair<double, double> p = std::make_pair(x0, y0);

	auto pos = std::find(points.begin(), points.end(), p);

	if (pos == points.end())
		return;

	this->N--;
	points.erase(pos);

	// recompute diffs
	this->diffs = divided_differences(this->points);
}

complex newton_poly::operator()(const complex &x) {
	complex val = complex(0.0, 0.0);

	for (size_t j = 0; j < N; j++) {
		val += diffs[j][0] * this->basis(x, j);
	}

	return val;
}

double newton_poly::operator()(const double x) {
	double val = 0;

	for (size_t j = 0; j < N; j++) {
		val += diffs[j][0] * this->basis(x, j);
	}

	return val;
}

// prints polynomial p in newton divided differences form.
std::ostream& operator<<(std::ostream& out, newton_poly &p) {
	for (size_t k = 0; k < p.N; k++) {
		double coeff = p.diffs[k][0];
		out << std::format("{}{}", (k == 0 ? "" : coeff >= 0 ? " + " : " - "), ABS(coeff));
		for (size_t j = 0; j < k; j++) {
			out << std::format("(x-{})", p.points[j].first);
		}
	}

	return out;
}
