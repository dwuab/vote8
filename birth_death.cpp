/*
 * birth_death.cpp
 *
 *  Created on: Nov 2, 2012
 *      Author: wudegang
 */

#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <cmath>
#include "boost/multi_array.hpp"

using namespace std;

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

double birth(int M) {

	double answer = 0.;
	double pt = p * eta / N;
	double qt = q * eta / N;
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
		for (int kp = 0; kp <= N - M; kp++)
			answer += d * combx(N - M, kp) * pow(qt, kp)
					* pow((1 - qt), N - M - kp) * heaviside(kp - k);
	}
	answer *= M;
	answer /= N;

	return answer;
}

double birth3(int x, int y, int z) {
	double answer = 0.;
	return answer;
}

double death3(int x, int y, int z) {
	double answer = 0.;
	return answer;
}

int main() {

	for (index i = 0; i < N + 1; ++i)
		for (index j = 0; j < N + 1; ++j)
			comb_table[i][j] = 0.;

	for (int i = 0; i < 901; i += 10)
	{
		cout << death(i) << ',';
	}
	cout << endl;
	return 0;
}
