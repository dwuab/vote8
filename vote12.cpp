// opinion formation on adaptive networks with intensive average degree with two communities (c-state Potts model)

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
#include <utility>
#include <list>

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

struct pairSecondGreater {
	bool operator()(const pair<int, int> *p1, const pair<int, int> *p2) {
		return (p1->second > p2->second);
	}
};

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << "opinion formation on adaptive networks with intensive average degree with c-state Potts model\n"
		     << "vote12 p eta\n";
		exit(1);
	}

	const int N = 900;
	const int Nrun = 1;
	const int c = 2;
	const double norm = 1. / (pow((c - 1.) / c, 2) + (c - 1) * pow(1. / c, 2));

	double p = boost::lexical_cast<double>(argv[1]);
	int eta = boost::lexical_cast<int>(argv[2]);

	double p_t = 1.0 * eta * p / N;
	double q = 1 - p;
	double q_t = 1.0 * eta * q / N;

	const int time_end = 25;
	const int time_point = 1 + time_end * 2;

	boost::random::uniform_int_distribution<> genN(0, N - 1);
	boost::random::uniform_int_distribution<> gen01(0, 1);
	boost::random::uniform_int_distribution<> gen3(0, 2);
	boost::random::uniform_real_distribution<> genReal;

	boost::mt19937 rng;
	rng.seed(static_cast<unsigned int>(std::time(0)));

	boost::array<int, N> spins;

	ofstream fout(
			str(format("N%d_p_%.3f_eta_%d_c_%d_T.dat") % N % p % c % eta).c_str());

	ofstream fout1(
			str(format("N%d_p_%.3f_eta_%d_c_%d_stat.dat") % N % p % c % eta).c_str());

	boost::array<int, c> M;
	bar randx(rng);

	list<pair<int, int>*> count_list;
	list<pair<int, int>*>::const_iterator cit = count_list.begin();
	list<pair<int, int>*>::iterator it = count_list.begin();
	for (int j = 0; j < c; ++j) {
		pair<int, int>* pair_ptr = new pair<int, int>(0,0);
		count_list.push_back(pair_ptr);
	}

	for (int i = 0; i < Nrun; ++i) {
		for (int j = 0; j < c; ++j)
			M[j] = 0;

		// initialization
		// in order to have exactly 1/c distribution, random shuffle is used
		for (int j = 0; j < c; ++j) {
			for (int k = (N / c) * j; k < (N / c) * (j + 1); ++k) {
				spins[k] = j;
			}
			M[j] = N / c;
		}

		random_shuffle(spins.begin(), spins.end(), randx);
//		for (int j = 0; j < N; ++j) {
//			spins[j] = gen3(rng);
//			++M[spins[j]];
//		}

		boost::array<double, c> x;
		for (int j = 0; j < c; ++j)
			x[j] = M[j] * 1.0 / N - 1.0 / c;

		for (int t = 0; t < N * time_end; ++t) {
			int pick = genN(rng);
			boost::array<int, c> N_count;
			for (int j = 0; j < c; ++j)
				N_count[j] = 0;

			for (int j = 0; j < N; ++j) {
				if (j == pick)
					continue;
				if (spins[j] == spins[pick])
					if (genReal(rng) < p_t) {
						N_count[spins[j]]++;

					} else {
					}
				else if (genReal(rng) < q_t) {
					N_count[spins[j]]++;

				}
			}

			// search for majority
			it = count_list.begin();
			for (int j=0;j<c;++j){
				(*it)->first=j;
				(*it)->second=N_count[j];
				++it;
			}

			count_list.sort(pairSecondGreater());
			cit = count_list.begin();
			int maj1, maj2;
			maj1 = (*cit)->second;
			++cit;
			maj2 = (*cit)->second;
			cit = count_list.begin();
			if (maj1 > maj2) {
				M[spins[pick]]--;
				M[(*cit)->first]++;
				spins[pick] = (*cit)->first;
			}

			// if there is no majority then we don't do anything

			if ((t + 1) % (N) == 0) {
//				for (int j = 0; j < c; ++j)
//					x[j] = M[j] * 1.0 / N - 1.0 / c;
//				for (int j = 0; j < c; ++j)
//					cout << M[j] << ' ';
//				cout << '\n';
			  int consensus = 0;
			  for (int j=0;j<c;++j)
			    consensus += static_cast<int>(M[j]==N);
			  if (consensus==1) //consensus reached
			    {
			      fout << (t+1)/N << '\n';
			      for (int j = 0; j < c; ++j)
				fout1 << M[j] << ' ';
			      break;
			    }
			  else if ((t+1)==N*time_end){
			    fout << time_end*10 << '\n';
			    for (int j = 0; j < c; ++j)
			      fout1 << M[j] << ' ';
			  }
			}

			if ((t + 1) == N * time_end) {
				for (int j = 0; j < c; ++j)
					x[j] = M[j] * 1.0 / N - 1.0 / c;
				double sq_sum = 0.;
				for (int j = 0; j < c; ++j)
					sq_sum += pow(x[j], 2);
				fout << sqrt(sq_sum * norm) << '\n';
				for (int j = 0; j < c; ++j)
					fout1 << M[j] << ' ';
				fout1 << '\n';

			}
		}
	}


	for (list<pair<int,int>* >::const_iterator it=count_list.begin();it!=count_list.end();++it)
		delete (*it);

	fout.close();
	fout1.close();

	return 0;
}
