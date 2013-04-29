/*
 * trans_prob.cpp
 *
 *  Created on: Apr 26, 2013
 *      Author: wudegang
 *
 *  Calculating the transition probabilities for 3-state model
 */

#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"
#include <boost/format.hpp>
#include <fstream>

using namespace std;
using boost::format;
using boost::io::str;

const int N = 150;
const double p = 0.5;
const double q = 1 - p;
const int eta = 10;

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index_t;

array_type comb_table(boost::extents[N + 1][N + 1]);
boost::array<double, N + 1> intermediate;

inline int heaviside(int x) {
	return static_cast<int>(x > 0);
}

inline double combx(int N, int k) {
	if (!(comb_table[N][k] > 0.)) {
		comb_table[N][k] = gsl_sf_choose(N, k);
	}
	return comb_table[N][k];
}

inline int strict_max(int x, int y) {
	if (x > y)
		return x;
	else if (y > x)
		return y;
	else
		return 0;
}

double W12(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int l = 0; l <= x - 1; ++l) {
		double c1 = combx(x - 1, l) * pow(pt, l) * pow(1 - pt, x - 1 - l);
		for (int lp = l + 1; lp <= y; ++lp) {
			double c2 = combx(y, lp) * pow(qt, lp) * pow(1 - qt, y - lp);
			for (int lpp = 0; lpp < lp && lpp <= z; ++lpp)
				intermediate[l] += c1 * c2 * combx(z, lpp) * pow(qt, lpp)
						* pow(1 - qt, z - lpp);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= x;
	result /= N;

	return result;
}

double W13(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int l = 0; l <= x - 1; ++l) {
		double c1 = combx(x - 1, l) * pow(pt, l) * pow(1 - pt, x - 1 - l);
		for (int lpp = l + 1; lpp <= z; ++lpp) {
			double c2 = combx(z, lpp) * pow(qt, lpp) * pow(1 - qt, z - lpp);
			for (int lp = 0; lp < lpp && lp <= y; ++lp)
				intermediate[l] += c1 * c2 * combx(y, lp) * pow(qt, lp)
						* pow(1 - qt, y - lp);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= x;
	result /= N;

	return result;
}

double W21(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int lp = 0; lp <= y - 1; ++lp) {
		double c1 = combx(y - 1, lp) * pow(pt, lp) * pow(1 - pt, y - 1 - lp);
		for (int l = lp + 1; l <= x; ++l) {
			double c2 = combx(x, l) * pow(qt, l) * pow(1 - qt, x - l);
			for (int lpp = 0; lpp < l && lpp <= z; ++lpp)
				intermediate[lp] += c1 * c2 * combx(z, lpp) * pow(qt, lpp)
						* pow(1 - qt, z - lpp);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= y;
	result /= N;

	return result;
}

double W23(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int lp = 0; lp <= y - 1; ++lp) {
		double c1 = combx(y - 1, lp) * pow(pt, lp) * pow(1 - pt, y - 1 - lp);
		for (int lpp = lp + 1; lpp <= z; ++lpp) {
			double c2 = combx(z, lpp) * pow(qt, lpp) * pow(1 - qt, z - lpp);
			for (int l = 0; l < lpp && l <= x; ++l)
				intermediate[lp] += c1 * c2 * combx(x, l) * pow(qt, l)
						* pow(1 - qt, x - l);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= y;
	result /= N;

	return result;
}

double W31(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int lpp = 0; lpp <= z - 1; ++lpp) {
		double c1 = combx(z - 1, lpp) * pow(pt, lpp) * pow(1 - pt, z - 1 - lpp);
		for (int l = lpp + 1; l <= x; ++l) {
			double c2 = combx(x, l) * pow(qt, l) * pow(1 - qt, x - l);
			for (int lp = 0; lp < l && lp <= y; ++lp)
				intermediate[lpp] += c1 * c2 * combx(y, lp) * pow(qt, lp)
						* pow(1 - qt, y - lp);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= z;
	result /= N;

	return result;
}

double W32(const int x, const int y) {
	double result = 0.0;
	double pt = p * eta / N;
	double qt = q * eta / N;
	const int z = N - x - y;

	for (int i = 0; i < N + 1; ++i) {
		intermediate[i] = 0.;
	}

#pragma omp parallel for
	for (int lpp = 0; lpp <= z - 1; ++lpp) {
		double c1 = combx(z - 1, lpp) * pow(pt, lpp) * pow(1 - pt, z - 1 - lpp);
		for (int lp = lpp + 1; lp <= y; ++lp) {
			double c2 = combx(y, lp) * pow(qt, lp) * pow(1 - qt, y - lp);
			for (int l = 0; l < lp && l <= x; ++l)
				intermediate[lpp] += c1 * c2 * combx(x, l) * pow(qt, l)
						* pow(1 - qt, x - l);
		}
	}

	for (int i = 0; i < N + 1; ++i) {
		result += intermediate[i];
	}

	result *= z;
	result /= N;

	return result;
}

int main() {
	ofstream fout12(
			str(
					format("data_trans_prob/N%d_p_%.3f_eta_%d_W12.dat") % N
							% p % eta).c_str());

	for (index_t i = 0; i < N + 1; ++i)
		for (index_t j = 0; j < N + 1; ++j)
			comb_table[i][j] = 0.;

	for (int i = 0; i < N + 1; i += 1) {
		cout << "i=" << i << endl;
		for (int j = 0; j <= N - i; j += 1) {

			fout12 << W12(i, j) << ' ';

		}
		fout12 << '\n';

	}

	fout12.close();

	return 0;
}
