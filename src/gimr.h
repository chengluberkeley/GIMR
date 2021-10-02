/**
 * Copyright 2016-2021, Cheng Lu, chenglu@berkeley.edu
 *
 * Header file of the algorithm:
 * "A faster algorithm solving a generalization of isotonic median regression
 * and a class of fused lasso problems", SIAM J. Optimization, 27(4), 2563-2596, 2017
 *
 * Author: Cheng Lu
 * Email: chenglu@berkeley.edu
 * Date: Oct. 9, 2016.
 */

#pragma once

#include "dp_array.h"
#include "dynamic_path.h"
#include <set>
#include <vector>

using namespace std;

struct ValueIndex {
	double value;
	int ext_index;  // Original index of the function 
	int int_index;  // Index of the breakpoint inside the function
	double left_slope;  // For integer solutions
	double right_slope;
};

struct SInterval {
	int left_node;
	int right_node;
};

// Comparison structure for the s_intervals
struct SIntervalComparison {
	bool operator() (const SInterval& s_interval1, const SInterval& s_interval2) const {
		return s_interval1.right_node < s_interval2.left_node;
	}
};

typedef set<SInterval, SIntervalComparison> s_interval_type;

/*
In our setup, to match the paper notation, we assume most arrays start from 1.
Thus the length of many arrays is 1 more than what's needed because of c++'s 0-offset convention.
*/
class GIMR {
  public:
	  // Public constructor
	  GIMR();

	  GIMR(int n, vector<vector<double>> w, vector<vector<double>> a, vector<double> d_ii1,
		   vector<double> d_i1i, bool is_integer);

	  int solve(clock_t t1, float diff);
	  int solve_dp();

	  void print();  // For the un-argument version, always apply for the internal data.
	  void print(const double* x);  // For outside arguments, always assume its range is from 0 to n-1.

	  double objective();
	  double objective(const double* x);

	  // Return the number of places of different values
	  int check_solution(const double* comp_x);

	  vector<double> x;  // Optimal solution

	  void reset();  // Reset the data structures ready to re-run the experiment.
  private:
	  // Function fields
	  int initialization_();
	  int initialization_dp_();
	  
	  // functions about the four arrays
	  int update_constant_(vector<double>& array, int i, double w);
	  int update_arrays_(int i_k, double w_ik_jk_1, double w_ik_jk);
	  int update_arrays_dp_(int i_k, double w_ik_jk_1, double w_ik_jk);
	  int find_status_change_interval_(int i_kl, int i_k, int i_kr, int& i_k1, int& i_k2);
	  int find_status_change_interval_dp_(int i_kl, int i_k, int i_kr, int& i_k1, int& i_k2);

	  // functions about s-intervals
	  int get_s_interval_(int i_k, int& i_kl, int& i_kr);
	  int update_s_interval_(int i_kl, int i_k1, int i_k2, int i_kr);

	  // Auxiliary function
	  double objective_(const double* x);
	  int locate_breakpoint_(double x, int index);
	  int positive_threshold_(double x);

	  // Function for integer cases
	  // Convert a function with fractional breakpoints into an equivalent one with integer breakpoints
	  int adjust_function_(const vector<double>& orig_w, const vector<double>& orig_a, const vector<double>& orig_b, vector<double>& new_w, vector<double>& new_a);

      // Data fields
	  // Input data fields
	  int n;  // Number of decision variables
	  vector<vector<double>> w;  // General weights
	  // vector<vector<double>> gurobi_w;
	  vector<vector<double>> a;  // General breakpoints
	  // vector<vector<double>> gurobi_a;
	  vector<double> d_ii1;  // Coefficients d_{i,i+1}
	  vector<double> d_i1i;  // Coefficients d_{i+1,i}
	  
	  // Other data structure
	  vector<ValueIndex> a_sorted;  // Sorted breakpoints
	  s_interval_type s_intervals;  // Red-black tree of s-intervals.
									// We use c++ stl "set" for red-black tree.
	  vector<double> sa, ta, tms, smt;  // Arrays that faciliates the computation 
									 // of node status changing intervals
	  vector<bool> status;  // False: status s, source set;
							// True: status t, sink set.
	  vector<vector<double>> b;  // The function values at the breakpoints.

	  // Dynamic paths
	  dp_array dp_sa, dp_ta, dp_tms, dp_smt;
	  dynamic_path_ops dp_ops;

	  // Marker for integer solution
	  bool is_integer;
};
