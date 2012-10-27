// opinion formation on adaptive networks with intensive average degree with two communities (3-state Potts model)

#define BOOST_DISABLE_ASSERTS // disable range-checking
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <functional>

using namespace std;
using boost::format;
using boost::io::str;

struct bar: std::unary_function<unsigned, unsigned> {
	boost::mt19937 &_state;
	unsigned operator()(unsigned i) {
		boost::uniform_int<> rng(0, i - 1);
		return rng(_state);
	}
	bar(boost::mt19937 &state) :
			_state(state) {
	}
};

int main(int argc, char** argv) {
	if (argc != 4) {
		cout
				<< "opinion formation on adaptive networks with intensive average degree with 3-state Potts model\n"
				<< "vote8 K p eta\n";
		exit(1);
	}

	const int N = 900;
	const int Nrun = 1;

	int K = boost::lexical_cast<int>(argv[1]);
	double p = boost::lexical_cast<double>(argv[2]);
	int eta = boost::lexical_cast<int>(argv[3]);

	double p_t_c = 1.0 * eta * p / N;
	double p_t = 1.0 * eta * p / K;
	double q = 1 - p;
	double q_t_c = 1.0 * eta * q / N;
	double q_t = 1.0 * eta * q / K;

	const int C_start = N - K;
	const int B_start = K;

	const int time_end = 500;
	const int time_point = 1 + time_end * 2;

	boost::random::uniform_int_distribution<> genN(0, N - 1);
	boost::random::uniform_int_distribution<> genK(0, K - 1);
	boost::random::uniform_int_distribution<> gen01(0, 1);
	boost::random::uniform_int_distribution<> gen3(0, 2);
	boost::random::uniform_real_distribution<> genReal;

	boost::mt19937 rng;
	rng.seed(static_cast<unsigned int>(std::time(0)));

	boost::array<int, N> spins;

	boost::multi_array<int, 2> adj_matrix(boost::extents[N][N]);

	ofstream fout(
			str(format("N%d_K_%d_p_%.3f_eta_%d_m.dat") % N % K % p % eta).c_str()
			);

	ofstream fout1(
			str(format("N%d_K_%d_p_%.3f_eta_%d_stat.dat") % N % K % p % eta).c_str()
			);

	boost::array<int, 3> M;
	boost::array<double, time_point> m;
	bar randx(rng);

	for (int i = 0; i < time_point; ++i)
		m[i] = 0.0;

	for (int i = 0; i < Nrun; ++i) {
		M[0] = 0;
		M[1] = 0;
		M[2] = 0;

		// initialization
		// in order to have exactly 1/3 1/3 1/3 distribution, random shuffle is used
		for (int j = 0; j < N / 3; ++j)
			spins[j] = 0;
		for (int j = N / 3; j < 2 * N / 3; ++j)
			spins[j] = 1;
		for (int j = 2 * N / 3; j < N; ++j)
			spins[j] = 2;
		M[0] = N / 3;
		M[1] = N / 3;
		M[2] = N / 3;

		random_shuffle(spins.begin(), spins.end(), randx);
//		for (int j = 0; j < N; ++j) {
//			spins[j] = gen3(rng);
//			++M[spins[j]];
//		}

		double x1, x2, x3;
		x1 = M[0] * 1.0 / N - 1.0 / 3;
		x2 = M[1] * 1.0 / N - 1.0 / 3;
		x3 = M[2] * 1.0 / N - 1.0 / 3;
//		m[0] += sqrt((pow(x1, 2) + pow(x2, 2) + pow(x3, 2)) / 3.) / Nrun;

		for (int t = 0; t < N * time_end; ++t) {
			int pick = genN(rng);
			boost::array<int, 3> N_count;
			N_count[0] = 0;
			N_count[1] = 0;
			N_count[2] = 0;
			int j_start = 0;
			int j_end = 0;
			double px = 0.0;
			double qx = 0.0;

			if (pick < C_start) // block A'
					{
//				cout << "A' ";
				j_start = 0;
				j_end = C_start;
				px = p_t;
				qx = q_t;

			} else if (pick >= C_start && pick < B_start) // block C
					{
//				cout << "C ";
				j_start = C_start;
				j_end = B_start;
				px = p_t_c;
				qx = q_t_c;
			} else { // block B'
//				cout << "B' ";
				j_start = B_start;
				j_end = N;
				px = p_t;
				qx = q_t;
			}

			for (int j = j_start; j < j_end; ++j) {
				if (j == pick)
					continue;
				if (spins[j] == spins[pick])
					if (genReal(rng) < px) {
						N_count[spins[j]]++;

					} else {
					}
				else if (genReal(rng) < qx) {
					N_count[spins[j]]++;

				}
			}

			// search for majority
			if (N_count[0] > N_count[1] && N_count[0] > N_count[2]) // spin 1 is majority
					{
				M[spins[pick]]--;
				M[0]++;
				spins[pick] = 0;
			} else if (N_count[1] > N_count[0] && N_count[1] > N_count[2]) // spin 2 is majority
					{
				M[spins[pick]]--;
				M[1]++;
				spins[pick] = 1;
			} else if (N_count[2] > N_count[0] && N_count[2] > N_count[1]) // spin 3 is majority
					{
				M[spins[pick]]--;
				M[2]++;
				spins[pick] = 2;
			}

			// if there is no majority then we don't do anything

			if ((t + 1) % (N) == 0) {
				x1 = M[0] * 1.0 / N - 1.0 / 3;
				x2 = M[1] * 1.0 / N - 1.0 / 3;
				x3 = M[2] * 1.0 / N - 1.0 / 3;
//				m[(t + 1) / (N / 2)] += sqrt(
//						(pow(x1, 2) + pow(x2, 2) + pow(x3, 2)) * 3.0 / 2.0)
//						/ Nrun;

//				cout << M[0] * 1.0 / N << ' ' << M[1] * 1.0 / N << ' '
//						<< M[2] * 1.0 / N << '\n';
//				cout << t << ' '
//						<< sqrt(
//								(pow(x1, 2) + pow(x2, 2) + pow(x3, 2)) * 3.0
//										/ 2.0) << '\n';
			}

			if ((t + 1) == N * time_end) {
				fout << sqrt((pow(x1, 2) + pow(x2, 2) + pow(x3, 2)) * 3.0 / 2.0)
						<< '\n';
				fout1 << M[0] << ' ' << M[1] << ' ' << M[2] << '\n';
//				for (int ii=0;ii<C_start;++ii)
//					cout << spins[ii] << ' ';
//				cout << '\n';
//				for (int ii=C_start;ii<B_start;++ii)
//					cout << spins[ii] << ' ';
//				cout << '\n';
//				for (int ii=B_start;ii<N;++ii)
//					cout << spins[ii] << ' ';
//				cout << '\n';
			}
		}
	}

//	for (int i = 0; i < time_point; i++)
//		cout << i * 0.5 << ' ' << m[i] << '\n';

	fout.close();
	fout1.close();

	return 0;
}
