/**
 * Copyright 2016-2021, Cheng Lu, chenglu@berkeley.edu

 * Implementation of the functions in gimr.h
 */

#include "gimr.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <utility>

// Comparison function for the breakpoints
bool breakpoint_comparison(const ValueIndex value_index1, const ValueIndex value_index2) {
	return value_index1.value < value_index2.value;
}

// Public constructor
GIMR::GIMR(int n, vector<vector<double>> w, vector<vector<double>> a,
	vector<double> d_ii1, vector<double> d_i1i, bool is_integer): 
	n(n), w(w), a(a), d_ii1(d_ii1), d_i1i(d_i1i), is_integer(is_integer) {
	// Check the validity of the data
	assert(n >= 1);

	// Sanity check of sizes
	assert(w.size() == (n + 1));
	assert(a.size() == (n + 1));
	assert(d_ii1.size() == (n + 1));
	assert(d_i1i.size() == (n + 1));

	for (int i = 1; i <= n; ++i) {
		assert(w[i].size() >= 2);  // Make sure at least two pieces
		assert(a[i].size() == w[i].size());  // pieces match breakpoints
		for (int j = 0; j < w[i].size() - 1; ++j)
			assert(w[i][j] < w[i][j + 1]);  // Convexity
		for (int j = 1; j < a[i].size() - 1; ++j)
			assert(a[i][j] < a[i][j + 1]);
	}

	// Check the validity of the separation coefficients
	assert((d_ii1[0] == 0) && (d_ii1[n] == 0) &&
		(d_i1i[0] == 0) && (d_i1i[n] == 0));
	for (int i = 1; i <= n - 1; ++i) {
		assert(d_ii1[i] >= 0);
		assert(d_i1i[i] >= 0);
	}

	// If all the above checks pass, initiate the x vector
	x.resize(n + 1, 0);

	// Compute the b vector for integer solution case
	// We assume objective value of the first breakpoint of 
	// each deviation function is 0. 
	// Solve the objective value of each breakpoint
	if (is_integer) {
		b.push_back(vector<double>(1, 0));
		for (int i = 1; i <= n; ++i) {
			vector<double> b_seq;
			b_seq.push_back(0);
			b_seq.push_back(0);
			for (int j = 2; j < a[i].size(); ++j) {
				double b_value = w[i][j - 1] * (a[i][j] - a[i][j - 1]) + b_seq[j - 1];
				b_seq.push_back(b_value);
			}
			b.push_back(b_seq);
		}
	}
}

GIMR::~GIMR() {
    clear_dps_();
}

// Public function members
int GIMR::solve(clock_t t1, float diff) {
	clock_t t2;
	initialization_();

	for (int k = 0; k < a_sorted.size(); ++k) {
		int i_k = a_sorted[k].ext_index;
		int j_k = a_sorted[k].int_index;		
		// Update graph
		if (!status[i_k]) {
			if (is_integer)
				update_arrays_(i_k, a_sorted[k].left_slope, a_sorted[k].right_slope);
			else
				update_arrays_(i_k, w[i_k][j_k - 1], w[i_k][j_k]);
		}

		if (!status[i_k]) {
			int i_kl, i_kr;
			get_s_interval_(i_k, i_kl, i_kr);
			int i_k1, i_k2;
			find_status_change_interval_(i_kl, i_k, i_kr, i_k1, i_k2);
			if (i_k1 > 0) {	
				for (int i = i_k1; i <= i_k2; ++i) {		
					x[i] = a_sorted[k].value;
					status[i] = true;
				}
				update_s_interval_(i_kl, i_k1, i_k2, i_kr);
			}
		}

		if (diff > 0) {
			t2 = clock();
			if (((float)t2 - (float)t1) / CLOCKS_PER_SEC > diff)
				return -1;
		}
	}

	return  0;
}

int GIMR::solve_dp() {
	initialization_dp_();
	
	for (int k = 0; k < a_sorted.size(); ++k) {
		int i_k = a_sorted[k].ext_index;
		int j_k = a_sorted[k].int_index;
		// Update graph
		if (!status[i_k]) {
			if (is_integer)
				update_arrays_dp_(i_k, a_sorted[k].left_slope, a_sorted[k].right_slope);
			else
				update_arrays_dp_(i_k, w[i_k][j_k - 1], w[i_k][j_k]);
		}

		if (!status[i_k]) {
			int i_kl, i_kr;
			get_s_interval_(i_k, i_kl, i_kr);
			int i_k1, i_k2;
			find_status_change_interval_dp_(i_kl, i_k, i_kr, i_k1, i_k2);
			if (i_k1 > 0) {
				for (int i = i_k1; i <= i_k2; ++i) {
					x[i] = a_sorted[k].value;
					status[i] = true;
				}
				update_s_interval_(i_kl, i_k1, i_k2, i_kr);
			}
		}
	}

	return  0;
}

void GIMR::print() {
	cout << "The optimal solution is: ";
	for (int i = 1; i <= n; ++i) {
		cout << x[i] << ' ';
	}
	cout << endl;
}

void GIMR::print(const double* x) {
	cout << "The optimal solution is: ";
	for (int i = 0; i < n; ++i) {
		cout << *(x + i) << ' ';
	}
	cout << endl;
}

double GIMR::objective() {
	double* conv_x = &x[0];
	return objective_(conv_x + 1);
}

double GIMR::objective(const double* x) {
	return objective_(x);
}

int GIMR::check_solution(const double* comp_x) {
	int value = 0;
	double max_diff = 0;
	for (int i = 1; i <= n; ++i) {
		// Special processing for double type data
		if (fabs(x[i] - *(comp_x + (i - 1))) >= 1e-6) {
			++value;
			if (fabs(x[i] - *(comp_x + (i - 1))) > max_diff)
				max_diff = fabs(x[i] - *(comp_x + (i - 1)));
		}
	}

	cout << "The maximum difference is " << max_diff << endl;

	return value;
}

void GIMR::reset() {
	a_sorted.clear();
	s_intervals.clear();
	sa.clear();
	ta.clear();
	tms.clear();
	smt.clear();
	status.clear();
    clear_dps_();
}

// Private function members
int GIMR::initialization_() {
	// Initialize the s_interval to be [1,n]
	SInterval g0_s_interval;
	g0_s_interval.left_node = 1;
	g0_s_interval.right_node = n;
	s_intervals.insert(g0_s_interval);

	// Initialize the four arrays
	sa.resize(n + 1, 0);
	ta.resize(n + 1, 0);
	tms.resize(n + 1, 0);
	smt.resize(n + 1, 0);

	for (int i = 1; i <= n; ++i) {
		sa[i] = sa[i - 1] - w[i][0];
	}

	for (int i = 1; i <= n; ++i) {
		tms[i] = ta[i] - sa[i] + d_ii1[i];
		smt[i] = sa[i] - ta[i] + d_i1i[i];
	}

	// Initialize the status array
	status.resize(n + 1, false);

	// Populate and sort the breakpoints
	if (is_integer) {
		ValueIndex value_index;
		for (int i = 1; i <= n; ++i) {
			int floor_a = floor(a[i][1]);
			int curr_bp = floor_a;
			double temp_b = w[i][0] * (floor_a - a[i][1]);
			value_index.value = floor_a;
			value_index.ext_index = i;
			value_index.left_slope = w[i][0];
			int ceil_a = ceil(a[i][1]);
			if (ceil_a > floor_a) {
				if ((w[i].size() == 2) || (ceil_a < a[i][2])) {
					double obj_value = w[i][1] * (ceil_a - a[i][1]);
					value_index.right_slope = obj_value - temp_b;
					a_sorted.push_back(value_index);
					value_index.value = ceil_a;
					value_index.left_slope = obj_value - temp_b;
					curr_bp = ceil_a;
					temp_b = obj_value;
				}
			}

			for (int j = 2; j < a[i].size(); ++j) {	
				floor_a = floor(a[i][j]);
				ceil_a = ceil(a[i][j]);
				if (ceil_a == floor_a) {
					value_index.right_slope = (b[i][j] - temp_b) / (a[i][j] - curr_bp);
					a_sorted.push_back(value_index);
					value_index.value = a[i][j];
					value_index.left_slope = value_index.right_slope;
					curr_bp = floor_a;
					temp_b = b[i][j];
				}
				else {
					if (floor_a > curr_bp) {
						value_index.right_slope = w[i][j - 1];
						a_sorted.push_back(value_index);
						value_index.value = floor_a;
						value_index.left_slope = value_index.right_slope;
						curr_bp = floor_a;
						temp_b = w[i][j - 1] * (floor_a - a[i][j]) + b[i][j];
					}

					if ((j == a[i].size() - 1) || (ceil_a < a[i][j + 1])) {
						double obj_value = w[i][j] * (ceil_a - a[i][j]) + b[i][j];
						value_index.right_slope = (obj_value - temp_b) / (ceil_a - curr_bp);
						a_sorted.push_back(value_index);
						value_index.value = ceil_a;
						value_index.left_slope = value_index.right_slope;
						curr_bp = ceil_a;
						temp_b = obj_value;
					}
				}
			}

			value_index.right_slope = w[i][w[i].size() - 1];
			a_sorted.push_back(value_index);
		}
	}
	else {
		for (int i = 1; i <= n; ++i)
			for (int j = 1; j < a[i].size(); ++j) {
				ValueIndex value_index;
				value_index.value = a[i][j];
				value_index.ext_index = i;
				value_index.int_index = j;
				a_sorted.push_back(value_index);
			}
	}

	sort(a_sorted.begin(), a_sorted.end(), breakpoint_comparison);

	return 0;
}

int GIMR::initialization_dp_() {
    // Clear existing ones.
    clear_dps_();

	// Initialize the s_interval to be [1,n]
	SInterval g0_s_interval;
	g0_s_interval.left_node = 1;
	g0_s_interval.right_node = n;
	s_intervals.insert(g0_s_interval);

	// Initialize the four arrays
	sa.resize(n + 1, 0);
	ta.resize(n + 1, 0);
	tms.resize(n + 1, 0);
	smt.resize(n + 1, 0);

	for (int i = 1; i <= n; ++i) {
		sa[i] = sa[i - 1] - w[i][0];
	}

    dp_sa = new dp_array(sa);
    dp_ta = new dp_array(ta);

	for (int i = 1; i <= n; ++i) {
		tms[i] = ta[i] - sa[i] + d_ii1[i];
		smt[i] = sa[i] - ta[i] + d_i1i[i];
	}

    dp_tms = new dp_array(tms);
    dp_smt = new dp_array(smt);

	// Initialize the status array
	status.resize(n + 1, false);

	// Populate and sort the breakpoints
	if (is_integer) {
		ValueIndex value_index;
		for (int i = 1; i <= n; ++i) {
			int floor_a = floor(a[i][1]);
			int curr_bp = floor_a;
			double temp_b = w[i][0] * (floor_a - a[i][1]);
			value_index.value = floor_a;
			value_index.ext_index = i;
			value_index.left_slope = w[i][0];
			int ceil_a = ceil(a[i][1]);
			if (ceil_a > floor_a) {
				if ((w[i].size() == 2) || (ceil_a < a[i][2])) {
					double obj_value = w[i][1] * (ceil_a - a[i][1]);
					value_index.right_slope = obj_value - temp_b;
					a_sorted.push_back(value_index);
					value_index.value = ceil_a;
					value_index.left_slope = obj_value - temp_b;
					curr_bp = ceil_a;
					temp_b = obj_value;
				}
			}

			for (int j = 2; j < a[i].size(); ++j) {
				floor_a = floor(a[i][j]);
				ceil_a = ceil(a[i][j]);
				if (ceil_a == floor_a) {
					value_index.right_slope = (b[i][j] - temp_b) / (a[i][j] - curr_bp);
					a_sorted.push_back(value_index);
					value_index.value = a[i][j];
					value_index.left_slope = value_index.right_slope;
					curr_bp = floor_a;
					temp_b = b[i][j];
				}
				else {
					if (floor_a > curr_bp) {
						value_index.right_slope = w[i][j - 1];
						a_sorted.push_back(value_index);
						value_index.value = floor_a;
						value_index.left_slope = value_index.right_slope;
						curr_bp = floor_a;
						temp_b = w[i][j - 1] * (floor_a - a[i][j]) + b[i][j];
					}

					if ((j == a[i].size() - 1) || (ceil_a < a[i][j + 1])) {
						double obj_value = w[i][j] * (ceil_a - a[i][j]) + b[i][j];
						value_index.right_slope = (obj_value - temp_b) / (ceil_a - curr_bp);
						a_sorted.push_back(value_index);
						value_index.value = ceil_a;
						value_index.left_slope = value_index.right_slope;
						curr_bp = ceil_a;
						temp_b = obj_value;
					}
				}
			}

			value_index.right_slope = w[i][w[i].size() - 1];
			a_sorted.push_back(value_index);
		}
	}
	else {
		for (int i = 1; i <= n; ++i)
			for (int j = 1; j < a[i].size(); ++j) {
				ValueIndex value_index;
				value_index.value = a[i][j];
				value_index.ext_index = i;
				value_index.int_index = j;
				a_sorted.push_back(value_index);
			}
	}

	sort(a_sorted.begin(), a_sorted.end(), breakpoint_comparison);

	return 0;
}

int GIMR::update_constant_(vector<double>& array, int i, double w) {
	for (int j = i; j < array.size(); ++j)
		array[j] = array[j] + w;
	return 0;
}

int GIMR::update_arrays_(int i_k, double w_ik_jk_1, double w_ik_jk) {
	if ((w_ik_jk_1 <= 0) && (w_ik_jk <= 0)) {
		update_constant_(sa, i_k, -(w_ik_jk - w_ik_jk_1));
		update_constant_(tms, i_k, w_ik_jk - w_ik_jk_1);
		update_constant_(smt, i_k, -(w_ik_jk - w_ik_jk_1));
	}
	else if ((w_ik_jk_1 <= 0) && (w_ik_jk >= 0)) {
		update_constant_(sa, i_k, w_ik_jk_1);
		update_constant_(ta, i_k, w_ik_jk);
		update_constant_(tms, i_k, w_ik_jk - w_ik_jk_1);
		update_constant_(smt, i_k, -(w_ik_jk - w_ik_jk_1));
	}
	else if ((w_ik_jk_1 >= 0) && (w_ik_jk >= 0)) {
		update_constant_(ta, i_k, w_ik_jk - w_ik_jk_1);
		update_constant_(tms, i_k, w_ik_jk - w_ik_jk_1);
		update_constant_(smt, i_k, -(w_ik_jk - w_ik_jk_1));
	}
	return 0;
}

int GIMR::update_arrays_dp_(int i_k, double w_ik_jk_1, double w_ik_jk) {
	if ((w_ik_jk_1 <= 0) && (w_ik_jk <= 0)) {
		dp_sa->update_constant(i_k, -(w_ik_jk - w_ik_jk_1));
		dp_tms->update_constant(i_k, w_ik_jk - w_ik_jk_1);
		dp_smt->update_constant(i_k, -(w_ik_jk - w_ik_jk_1));
	}
	else if ((w_ik_jk_1 <= 0) && (w_ik_jk >= 0)) {
		dp_sa->update_constant(i_k, w_ik_jk_1);
		dp_ta->update_constant(i_k, w_ik_jk);
		dp_tms->update_constant(i_k, w_ik_jk - w_ik_jk_1);
		dp_smt->update_constant(i_k, -(w_ik_jk - w_ik_jk_1));
	}
	else if ((w_ik_jk_1 >= 0) && (w_ik_jk >= 0)) {
		dp_ta->update_constant(i_k, w_ik_jk - w_ik_jk_1);
		dp_tms->update_constant(i_k, w_ik_jk - w_ik_jk_1);
		dp_smt->update_constant(i_k, -(w_ik_jk - w_ik_jk_1));
	}
	return 0;
}

int GIMR::find_status_change_interval_(int i_kl, int i_k, int i_kr, int& i_k1, int& i_k2) {
	// Find the optimal solution when the status change interval is not empty
	int hat_i_k1 = i_kl;
	double f1_hat = sa[i_k] - sa[i_kl - 1];
	if (i_kl + 1 <= i_k) {
		// Find tilde{i_{k1}}
		int tilde_i_k1 = i_kl + 1;
		double curr_min = tms[i_kl];
		for (int i = i_kl + 2; i <= i_k; ++i) {
			if (tms[i - 1] <= curr_min) {
				curr_min = tms[i - 1];
				tilde_i_k1 = i;
			}
		}

		// Find hat{i_{k1}}
		double f1_tilde_i_k1 = tms[tilde_i_k1 - 1] - ta[i_kl - 1] + sa[i_k] + d_i1i[i_kl - 1];

		if (f1_tilde_i_k1 <= f1_hat) {
			hat_i_k1 = tilde_i_k1;
			f1_hat = f1_tilde_i_k1;
		}
	}

	int hat_i_k2 = i_kr;
	double f2_hat = sa[i_kr] - sa[i_k];
	if (i_k <= i_kr - 1) {
		// Find tilde{i_{k2}}
		int tilde_i_k2 = i_kr - 1;
		double curr_min = smt[i_kr - 1];
		for (int i = i_kr - 2; i >= i_k; --i) {
			if (smt[i] <= curr_min) {
				curr_min = smt[i];
				tilde_i_k2 = i;
			}
		}

		// Find hat{i_{k2}}
		double f2_tilde_i_k2 = smt[tilde_i_k2] - sa[i_k] + ta[i_kr] + d_ii1[i_kr];

		if (f2_tilde_i_k2 <= f2_hat) {
			hat_i_k2 = tilde_i_k2;
			f2_hat = f2_tilde_i_k2;
		}
	}

	// Compute the objective of empty status change interval
	double z_empty = ta[i_kr] - ta[i_kl - 1] + d_i1i[i_kl - 1] + d_ii1[i_kr];

	if (f1_hat + f2_hat >= z_empty) {
		i_k1 = 0;
		i_k2 = 0;
	}
	else {
		i_k1 = hat_i_k1;
		i_k2 = hat_i_k2;
	}

	return 0;
}

int GIMR::find_status_change_interval_dp_(int i_kl, int i_k, int i_kr, int& i_k1, int& i_k2) {
	// Precompute some terms
    double sa_i_k = dp_sa->edge_cost(i_k);
    double sa_i_kl_1 = dp_sa->edge_cost(i_kl - 1);
    double sa_i_kr = dp_sa->edge_cost(i_kr);
    double ta_i_kl_1 = dp_ta->edge_cost(i_kl - 1);
    double ta_i_kr = dp_ta->edge_cost(i_kr);
	// Find the optimal solution when the status change interval is not empty
	int hat_i_k1 = i_kl;
	double f1_hat =  sa_i_k - sa_i_kl_1;
	if (i_kl + 1 <= i_k) {
		// Find tilde{i_{k1}}
		int tilde_i_k1;
        double cost = dp_tms->min_cost_last(i_kl, i_k, tilde_i_k1);

		// Find hat{i_{k1}}
		double f1_tilde_i_k1 = cost - ta_i_kl_1 + sa_i_k + d_i1i[i_kl - 1];

		if (f1_tilde_i_k1 <= f1_hat) {
			hat_i_k1 = tilde_i_k1;
			f1_hat = f1_tilde_i_k1;
		}
	}

	int hat_i_k2 = i_kr;
	double f2_hat = sa_i_kr - sa_i_k;
	if (i_k <= i_kr - 1) {
		// Find tilde{i_{k2}}
		int tilde_i_k2;
        double cost = dp_smt->min_cost_first(i_k, i_kr, tilde_i_k2);

		// Find hat{i_{k2}}
		double f2_tilde_i_k2 = cost - sa_i_k + ta_i_kr + d_ii1[i_kr];

		if (f2_tilde_i_k2 <= f2_hat) {
			hat_i_k2 = tilde_i_k2;
			f2_hat = f2_tilde_i_k2;
		}
	}

	// Compute the objective of empty status change interval
	double z_empty = ta_i_kr - ta_i_kl_1 + d_i1i[i_kl - 1] + d_ii1[i_kr];

	if (z_empty <= f1_hat + f2_hat) {
		i_k1 = 0;
		i_k2 = 0;
	}
	else {
		i_k1 = hat_i_k1;
		i_k2 = hat_i_k2;
	}

	return 0;
}

int GIMR::get_s_interval_(int i_k, int& i_kl, int& i_kr) {
	SInterval temp_s_interval;
	temp_s_interval.left_node = i_k;
	temp_s_interval.right_node = i_k;
	s_interval_type::iterator it = s_intervals.upper_bound(temp_s_interval);
	--it;  // Get the actual iterator containing the node i_k
	i_kl = it->left_node;
	i_kr = it->right_node;

	return 0;
}

int GIMR::update_s_interval_(int i_kl, int i_k1, int i_k2, int i_kr) {
	SInterval temp_s_interval;
	temp_s_interval.left_node = i_kl;
	temp_s_interval.right_node = i_kr;
	s_intervals.erase(temp_s_interval);

	if (i_kl <= i_k1 - 1) {
		temp_s_interval.left_node = i_kl;
		temp_s_interval.right_node = i_k1 - 1;
		s_intervals.insert(temp_s_interval);
	}

	if (i_k2 + 1 <= i_kr) {
		temp_s_interval.left_node = i_k2 + 1;
		temp_s_interval.right_node = i_kr;
		s_intervals.insert(temp_s_interval);
	}

	return 0;
}

// Assume the data starts at offset 0.
double GIMR::objective_(const double* x) {
	double obj = 0;
	for (int i = 1; i <= n; ++i) {
		int loc = locate_breakpoint_(*(x+(i-1)), i);
		double y_value;
		if (loc == 0) {
			y_value = w[i][0] * (*(x+(i-1)) - a[i][1]);
			obj += y_value;
		}
		else {
			y_value = w[i][loc] * (*(x+(i-1)) - a[i][loc]) + b[i][loc];
			obj += y_value;
		}
	}

	for (int i = 1; i <= n - 1; ++i) {
		obj += d_ii1[i] * positive_threshold_(*(x + (i - 1)) - *(x + i))
			+ d_i1i[i] * positive_threshold_(*(x + i) - *(x + (i - 1)));
	}

	return obj;
}

int GIMR::locate_breakpoint_(double x, int index) {
	if (x < a[index][1]) return 0;
	for (int i = 1; i < a[index].size() - 1; ++i) {
		if ((x >= a[index][i]) && (x < a[index][i + 1]))
			return i;
	}

	return a[index].size() - 1;
}

int GIMR::positive_threshold_(double x) {
	if (x > 0)
		return x;
	else
		return 0;
}

int GIMR::adjust_function_(const vector<double>& orig_w, const vector<double>& orig_a, const vector<double>& orig_b, vector<double>& w, vector<double>& a) {
	vector<double> temp_b;
	temp_b.push_back(0);
	w.clear();
	a.clear();
	a.push_back(0);  // Not used
	// Dealing with the first breakpoint
	if (orig_a[1] - floor(orig_a[1]) < 1e-6) { // Is integer
		a.push_back(orig_a[1]);
		w.push_back(orig_w[0]);
		temp_b.push_back(0);  // Assume the first breakpoint has function value 0.
	}
	else {
		int floor_a = floor(orig_a[1]);
		int ceil_a = ceil(orig_a[1]);
		a.push_back(floor_a);
		w.push_back(orig_w[0]);
		temp_b.push_back(orig_w[0] * (floor_a - orig_a[1])); // Assume the first breakpoint has function value 0.
		if ((orig_w.size() == 2) || (ceil_a < orig_a[2])) {
			// Insert ceil_a
			double obj_value = orig_w[1] * (ceil_a - orig_a[1]);
			a.push_back(ceil_a);
			w.push_back(obj_value - temp_b[1]);
			temp_b.push_back(obj_value);
		}
	}

	for (int i = 2; i < orig_a.size(); ++i) {
		int floor_a = floor(orig_a[i]);
		int ceil_a = ceil(orig_a[i]);
		if ((orig_a[i] - floor_a) < 1e-6) {  // Integer case
			w.push_back((orig_b[i] - temp_b[temp_b.size() - 1]) / (orig_a[i] - a[a.size() - 1]));
			a.push_back(orig_a[i]);
			temp_b.push_back(orig_b[i]);
		} else {
			if (floor_a > a[a.size() - 1]) {
				a.push_back(floor_a);
				w.push_back(orig_w[i - 1]);
				temp_b.push_back(orig_w[i - 1] * (floor_a - orig_a[i]) + orig_b[i]);
			}

			if ((i == orig_a.size() - 1) || (ceil_a < orig_a[i + 1])) {
				// Insert ceil_a
				double obj_value = orig_w[i] * (ceil_a - orig_a[i]) + orig_b[i];
				w.push_back((obj_value - temp_b[temp_b.size() - 1]) / (ceil_a - a[a.size() - 1]));
				a.push_back(ceil_a);
				temp_b.push_back(obj_value);
			}
		}	
	}

	// Insert the last slope
	w.push_back(orig_w[orig_w.size() - 1]);

	return 0;
}

void GIMR::clear_dps_() {
    if (dp_sa) {
        delete dp_sa;
        dp_sa = nullptr;
    }

    if (dp_ta) {
        delete dp_ta;
        dp_ta = nullptr;
    }

    if (dp_tms) {
        delete dp_tms;
        dp_tms = nullptr;
    }

    if (dp_smt) {
        delete dp_smt;
        dp_smt = nullptr;
    }
}
