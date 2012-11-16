// opinion formation on adaptive networks with intensive average degree with two communities

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
	if (argc != 3) {
		cout
				<< "opinion formation on adaptive networks with intensive average degree with community structure\n"
				<< "vote7 p eta\n";
		exit(1);
	}

	const int N = 900;
	const int Nrun = 10000;

	double p = boost::lexical_cast<double>(argv[1]);
	int eta = boost::lexical_cast<int>(argv[2]);

	const int M0 = 450;
	const double up_frac = 1.0 * M0 / N;
	const double p_t = 1.0 * eta * p / N;
	const double q = 1 - p;
	const double q_t = 1.0 * eta * q / N;

	const int time_end = 20000;
	const int time_point = 1 + time_end * 2;

	boost::random::uniform_int_distribution<> genN(0, N - 1);
	boost::random::uniform_int_distribution<> gen01(0, 1);
	boost::random::uniform_int_distribution<> gen3(0, 2);
	boost::random::uniform_real_distribution<> genReal;

	boost::mt19937 rng;
	rng.seed(static_cast<unsigned int>(std::time(0)));

	boost::array<int, N> spins;

	boost::multi_array<int, 2> adj_matrix(boost::extents[N][N]);

	ofstream fout(str(format("N%d_p_%.3f_eta_%d_T.dat") % N % p % eta).c_str());

	// this file records the opinion distribution at the end of the simulation
	ofstream fout1(
			str(format("N%d_p_%.3f_eta_%d_stat.dat") % N % p % eta).c_str());

	bar randx(rng);

	int N_tie = 0;

	for (int i = 0; i < Nrun; ++i) { // each loop is an independent simulation
		cout << i << endl;
		int M = 0;

		for (int j = 0; j < N / 2; ++j)
			spins[j] = 1;
		for (int j = N / 2; j < N; ++j)
			spins[j] = -1;
		M = N / 2;

//		m[0] += 1.0 * M / N / Nrun;

		for (int t = 0; t < N * time_end; ++t) {
			int pick = genN(rng);
			int N_nei = 0;

			int pCount = 0; // number of plus
			int mCount = 0; // number of minus

			for (int j = 0; j < N; ++j) {
				if (j == pick)
					continue;
				if (spins[j] == spins[pick])
					if (genReal(rng) < p_t) {
						if (spins[j] == 1)
							pCount++;
						else
							mCount++;
					} else {
					}
				else if (genReal(rng) < q_t) {
					if (spins[j] == 1)
						pCount++;
					else
						mCount++;
				}
			}

			// search for majority
			int spin_sum = pCount - mCount;
			if (spin_sum > 0) {
				if (spins[pick] == -1)
					++M;
				spins[pick] = 1;
			} else if (spin_sum < 0) {
				if (spins[pick] == 1)
					--M;
				spins[pick] = -1;
			} else
			//N_tie += 1;

			if ((M == 0) || (M == N)) {
				fout << 1.0 * (t + 1) / N << endl;
				fout1 << M << '\n';
				break;
			}

			if ((t + 1) % (N) == 0) {
				//fout1 << N_tie << ' ';
				//N_tie = 0;
//				m[(t + 1) / (N)] += 1.0 * M / N / Nrun;
//				cout << 1.0 * M / N << '\n';

			}
			if (t + 1 == N * time_end) {
				fout << 1.0 * (t + 1) / N << '\n';
				fout1 << M << '\n';
			}
		}
	}

	fout.close();
	fout1.close();

	return 0;
}
