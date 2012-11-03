/*
 * birth_death.cpp
 *
 *  Created on: Nov 2, 2012
 *      Author: wudegang
 */

#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "boost/multi_array.hpp"
#include <boost/format.hpp>
#include <fstream>

using namespace std;
using boost::format;
using boost::io::str;

const int N = 900;
const double p = 0.3;
const double q = 1 - p;
const int eta = 10;

typedef boost::multi_array<double, 2> array_type;
typedef array_type::index index;

array_type comb_table(boost::extents[N + 1][N + 1]);

inline int heaviside(int x) {
	return static_cast<int>(x > 0);
}

inline double combx(int N, int k) {
	if (!(comb_table[N][k] > 0.))
		comb_table[N][k] = gsl_sf_choose(N, k);
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

double birth(int M) {

	double answer = 0.;
	double pt = p * eta / N;
	double qt = q * eta / N;

#pragma omp parallel for
	for (int l = 0; l <= N - M - 1; ++l) {
		double c = combx(N - M - 1, l);
		double d = c * pow(pt, l) * pow((1 - pt), N - M - 1 - l);
		for (int lp = 0; lp <= M; lp++)
			answer += d * combx(M, lp) * pow(qt, lp) * pow((1 - qt), M - lp)
					* heaviside(lp - l);
	}
	answer *= N - M;
	answer /= N;

	return answer;
}

double death(int M) {
	double answer = 0.;
	double pt = p * eta / N;
	double qt = q * eta / N;
	for (int k = 0; k <= M - 1; ++k) {
		double c = combx(M - 1, k);
		double d = c * pow(pt, k) * pow((1 - pt), M - 1 - k);
		cout << d << endl;
		for (int kp = 0; kp <= N - M; kp++)
			if (kp > k)
				answer += d * combx(N - M, kp) * pow(qt, kp)
						* pow((1 - qt), N - M - kp);
	}
	answer *= M;
	answer /= N;

	return answer;
}

// birth rate calculation for 3 state model
// in the following comments, we assume there are x spin 1, y spin 2
// z (=N-x-y) spin 3
double birth3(const int x, const int y) {
	double answer1 = 0.;
	double answer2 = 0.;
	const double pt = p * eta / N;
	const double qt = q * eta / N;
	const int z = N - x - y;

	// suppose a spin 2 is chosen
#pragma omp parallel for
	for (int l = 0; l <= y - 1; ++l) {
		const double c1 = combx(y - 1, l) * pow(pt, l) * pow(1 - pt, y - 1 - l);
		for (int lpp = 0; lpp <= z; ++lpp) {
			const double c2 = combx(z, lpp) * pow(qt, lpp)
					* pow(1 - qt, z - lpp);
			for (int lp = 0; lp <= x; ++lp)
				if ((lp > l) && (lp > lpp))
					answer1 += c1 * c2 * combx(x, lp) * pow(qt, lp)
							* pow(1 - qt, x - lp);
		}
	}

	answer1 *= 1.0 * y / N;

	// suppose a spin 3 is chosen
	for (int l = 0; l <= z - 1; ++l) {
		const double c1 = combx(z - 1, l) * pow(pt, l) * pow(1 - pt, z - 1 - l);
		for (int lpp = 0; lpp <= y; ++lpp) {
			const double c2 = combx(y, lpp) * pow(qt, lpp)
					* pow(1 - qt, y - lpp);
			for (int lp = 0; lp <= x; ++lp)
				if ((lp > l) && (lp > lpp))
					answer2 += c1 * c2 * combx(x, lp) * pow(qt, lp)
							* pow(1 - qt, x - lp);
		}
	}

	answer2 *= 1.0 * z / N;

	return answer1 + answer2;
}

// death rate calculation for 3 state model
// in the following comments, we assume there are x spin 1, y spin 2
// z (=N-x-y) spin 3
double death3(int x, int y) {
	double answer = 0.;
	double pt = p * eta / N;
	double qt = q * eta / N;
	int z = N - x - y;

#pragma omp parallel for
	for (int k = 0; k <= x - 1; ++k) {
		double c1 = combx(x - 1, k) * pow(pt, k) * pow((1 - pt), x - 1 - k);
		for (int kp = 0; kp <= y; ++kp) {
			double c2 = combx(y, kp) * pow(qt, kp) * pow((1 - qt), y - kp);
			for (int kpp = 0; kpp <= z; ++kpp) {
				if (strict_max(kp, kpp) > k)
					answer += c1 * c2 * combx(z, kpp) * pow(qt, kpp)
							* pow((1 - qt), z - kpp);
			}
		}
	}

	answer *= x;
	answer /= N;

	return answer;
}

int main() {
	ofstream fout(
			str(format("N%d_p_%.3f_eta_%d_birth.dat") % N % p % eta).c_str());

	ofstream fout1(
			str(format("N%d_p_%.3f_eta_%d_death.dat") % N % p % eta).c_str());

	for (index i = 0; i < N + 1; ++i)
		for (index j = 0; j < N + 1; ++j)
			comb_table[i][j] = 0.;

	for (int i = 0; i < 901; i += 10) {
		cout << "i=" << i << endl;
		for (int j = 0; j <= N - i; j += 10) {
			fout << birth3(i, j) << ' ';
			fout1 << death3(i, j) << ' ';
		}
		fout << '\n';
		fout1 << '\n';
	}

	fout.close();
	fout1.close();
	return 0;
}
