/**
 * This code is an implementation of the algorithm of the following paper:
 * "A faster algorithm solving a generalization of isotonic median regression
 *  and a class of fused lasso problems", SIAM J. Optimization, 27(4), 2563-2596, 2017
 * We solve the following problem:
 * (GIMR) min_{x_1,...,x_n} \sum_{i=1}^n f^{pl}_i(x_i; \{a_{i,j}\}^{q_i}_{j=1})
 *	+ \sum_{i=1}^{n-1}d_{i,i+1}(x_i - x_{i+1})_+
 *	+ \sum_{i=1}^{n-1}d_{i+1,i}(x_{i+1} - x_i)_+
 * s.t. \ell_i \leq x_i \leq u_i, i = 1,...,n
 *
 * This is a sample main file to illustrate how to call our solver.
 *
 * Author: Cheng Lu
 * Email: chenglu@berkeley.edu
 * Date: Oct. 9, 2016.
 */

#include "gimr.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>

const int CONST_Q = 10000;
const int MODULA = 100;
const int MAX_INT = 5;
const bool IS_INTEGER = false;
const bool IS_BP_INTEGER = false;
const string FILENAME = "/tmp/result_10000.txt";


/*
We change everything to double to test the performance of 
different algorithms in finding integer solutions.
*/

vector<double> generate_random_sequence(int n, bool is_integer) {
	vector<double> seq(n, 0);
	for (int i = 0; i < n; ++i) {
		double value = rand() % MODULA + 1;
		if (is_integer)
			seq[i] = value;
		else
			seq[i] = value / (rand() % MODULA + 5);
	}

	return seq;
}

bool same_value(double a, double b) {
	double exp = 1e-6;
	if (std::fabs(a - b) < exp)
		return true;
	else
		return false;
}

bool duplicate(const vector<double> seq, double value, int index) {
	for (int i = 0; i < index; ++i) {
		if (same_value(seq[i], value))
			return true;
	}
	if (same_value(seq[seq.size() - 1], value))
		return true;

	return false;
}

vector<double> generate_random_sequence_no_repeat(int n, bool is_integer) {
	/* To facilitate the set-up of our algorithm, we shall guarantee 
	   that the first element of the sequence is negative and the last 
	   element of the sequence is positive. */

	vector<double> seq(n, 0);

	if (is_integer) {
		seq[0] = rand() % MODULA - MODULA;
		seq[n - 1] = rand() % MODULA + 1;
		// The middle element can be arbitrary sign.
		for (int i = 1; i < n - 1; ++i) {
			double value = rand() % MODULA - MODULA / 2;
			while (duplicate(seq, value, i))
				value = rand() % MODULA - MODULA / 2;
			seq[i] = value;
		}
	}
	else {
		seq[0] = ((double)(rand() % MODULA - MODULA)) / (rand() % MODULA + 5);
		seq[n - 1] = ((double)(rand() % MODULA + 1)) / (rand() % MODULA + 5);
		// The middle element can be arbitrary sign.
		for (int i = 1; i < n - 1; ++i) {
			double value = rand() % MODULA - MODULA / 2;
			value = value / (rand() % MODULA + 5);
			while (duplicate(seq, value, i))
				value = ((double)(rand() % MODULA - MODULA / 2)) / (rand() % MODULA + 5);
			seq[i] = value;
		}
	}

	return seq;
}

int generate_abs_data(int n, vector<int>& q, vector<vector<double>>& w,
	vector<vector<double>>& a, vector<double>& d_ii1,
	vector<double>& d_i1i) {
	// Clear the vectors for re-use
	q.clear();
	w.clear();
	a.clear();
	d_ii1.clear();
	d_i1i.clear();
	// Generate input data of weighted absolute deviation functions
	// srand(time(NULL));
	q.resize(n + 1, 1);
	w.push_back(vector<double>(1, 0));
	a.push_back(vector<double>(1, 0));
	for (int i = 1; i <= n; ++i) {
		vector<double> w_seq;
		vector<double> a_seq;
		vector<double> value = generate_random_sequence(1, IS_BP_INTEGER);
		w_seq.push_back(-value[0]);
		w_seq.push_back(value[0]);
		w.push_back(w_seq);
		a_seq.push_back(0);
		value = generate_random_sequence(1, IS_BP_INTEGER);
		a_seq.push_back(value[0]);
		a.push_back(a_seq);
	}
	d_ii1 = generate_random_sequence(n + 1, IS_BP_INTEGER);
	// d_i1i = generate_random_sequence(n + 1);
	d_i1i = d_ii1;
	d_ii1[0] = 0;
	d_ii1[n] = 0;
	d_i1i[0] = 0;
	d_i1i[n] = 0;

	cout << "Finish generating absolute value data." << endl;

	return 0;
}

int generate_data(int n, int bar_q, vector<int>& q, vector<vector<double>>& w,
	vector<vector<double>>& a, vector<double>& d_ii1,
	vector<double>& d_i1i) {
	// Clear the vectors for re-use
	q.clear();
	w.clear();
	a.clear();
	d_ii1.clear();
	d_i1i.clear();
	// Generate random input values
	// srand(time(NULL));
	// Generate the input data
	q.resize(n + 1, bar_q);  // We make it deterministic, the most extreme case 
	// int total_q = 0;
	// for (vector<int>::iterator it = q.begin()+1; it != q.end(); ++it)
	// 	total_q += *it;
	cout << "The total number of breakpoints is " << n*bar_q << endl;
	w.push_back(vector<double>(1, 0));
	a.push_back(vector<double>(1, 0));
	for (int i = 1; i <= n; ++i) {
		vector<double> w_seq = generate_random_sequence_no_repeat(q[i] + 1, IS_BP_INTEGER);
		vector<double> a_seq = generate_random_sequence_no_repeat(q[i] + 1, IS_BP_INTEGER);
		sort(w_seq.begin(), w_seq.end());
		sort(a_seq.begin(), a_seq.end());
		w.push_back(w_seq);
		a.push_back(a_seq);
	}
	d_ii1 = generate_random_sequence(n + 1, IS_BP_INTEGER);
	//d_i1i = generate_random_sequence(n + 1);
	d_i1i = d_ii1;
	d_ii1[0] = 0;
	d_ii1[n] = 0;
	d_i1i[0] = 0;
	d_i1i[n] = 0;

	cout << "Finish generating general piecewise linear data." << endl;

	return 0;
}

float mean(const vector<float>& v) {
	float value = 0;
	for (int i = 0; i < v.size(); ++i)
		value = value + v[i];
	if (v.size() > 0)
		return value / v.size();
	else
		return 0;
}

float variance(const vector<float>& v, float mean_value) {
	float value = 0;
	for (int i = 0; i < v.size(); ++i)
		value = value + (v[i] - mean_value) * (v[i] - mean_value);
	if (v.size() > 0)
		return value / v.size();
	else
		return 0;
}


int main() {
	// Try different pairs of (n, bar{q}). 
	// For each pair, 5 runs with random data
	int n, bar_q;
	bar_q = 1000;

	// Initialize the timers
	clock_t t1 = 0, t2 = 0;
	// Initialize file
	ofstream myfile;
	myfile.open(FILENAME, ios::out | ios::app);
	srand(time(NULL));
	while (bar_q < CONST_Q) {
		n = CONST_Q / bar_q;
		vector<float> our_diff_dp(MAX_INT, 0);
		vector<float> our_diff(MAX_INT, 0);

		// Generate data
		vector<int> q;
		vector<vector<double>> w, a;
		vector<double> d_ii1, d_i1i;  // d_{i,i+1}, d_{i+1,i}
		bool is_integer = IS_INTEGER;
		int num_early_stop = 0;

		for (int i = 0; i < MAX_INT; ++i) {
			if (bar_q == 1)
				generate_abs_data(n, q, w, a, d_ii1, d_i1i);
			else
				generate_data(n, bar_q, q, w, a, d_ii1, d_i1i);

			// Method 1: Our algorithm with dynamic path implementation	
			cout << "Round 1: Our algorithm (Dynamic paths)" << endl;
			t1 = clock();
			GIMR solver_dp(n, w, a, d_ii1, d_i1i, is_integer);
			solver_dp.solve_dp();
			t2 = clock();
			our_diff_dp[i] = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;
			cout << "The total running time of our algorithm (dynamic path) is " << our_diff_dp[i] << "s" << endl;
			cout << endl << endl;

			// Method 2: Our algorithm without dynamic path implementation	
			cout << "Round 2: Our algorithm (NO Dynamic paths)" << endl;
			int status;
			t1 = clock();
			GIMR solver(n, w, a, d_ii1, d_i1i, is_integer);
			if (n <= 10000)
				status = solver.solve(t1, -1); // Finish all computation
			else
				status = solver.solve(t1, our_diff_dp[i]); // Early stop
			t2 = clock();
			our_diff[i] = ((float)t2 - (float)t1) / CLOCKS_PER_SEC;

			if (status < 0) {
				num_early_stop++;
				cout << "Early stops!" << endl;
			}
			cout << "The total running time of our algorithm (NO dynamic path) is " << our_diff[i] << "s" << endl;
			cout << endl << endl;
		}

		// Compute the mean and variance
		float our_mean_dp, our_variance_dp, our_mean, our_variance;
		our_mean_dp = mean(our_diff_dp);
		our_variance_dp = variance(our_diff_dp, our_mean_dp);
		our_mean = mean(our_diff);
		our_variance = variance(our_diff, our_mean);

		cout << "The final aggregate results are as follows: " << endl;
		cout << "Our algorithm (dynamic path): mean = " << our_mean_dp << "s, variance = " << our_variance_dp << "s" << endl;
		cout << "Our algorithm (NO dynamic path): mean = " << our_mean << "s, variance = " << our_variance << "s" << endl;
		cout << "The number of early stops for NO dynamic path is " << num_early_stop << endl;
		cout << endl << endl;

		// Output one result to file
		myfile << "(n, bar_q) = (" << n << ", " << bar_q << ")" << endl;
		myfile << "The final aggregate results are as follows: " << endl;
		myfile << "Our algorithm (dynamic path): mean = " << our_mean_dp << "s, variance = " << our_variance_dp << "s" << endl;
		myfile << "Our algorithm (NO dynamic path): mean = " << our_mean << "s, variance = " << our_variance << "s" << endl;
		myfile << "The number of early stops for NO dynamic path is " << num_early_stop << endl;
		myfile << endl << endl;

		bar_q = bar_q * 10;
	}

	myfile.close();

	return 0;
}
