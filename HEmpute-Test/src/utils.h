/*
 * @file       utils.h, header file
 * @brief      defining functions for reading and storing data
 *
 * @author     Miran Kim
 * @date       June. 3, 2019
 * @copyright  GNU Pub License
 */

#ifndef UTILS_H_
#define UTILS_H_
 
#include <iostream>
#include <vector>
#include <string>

using namespace std;

typedef vector< vector< vector< vector<double> > > > ften;
typedef vector< vector< vector<double> > > dten;
typedef vector< vector<double> >  dmat;
typedef vector<double>  dvec;
typedef vector<uint64_t>  uvec;


bool scale_and_modt(uint64_t& x, double a, double scale, size_t coeff_modulus);
bool scale_and_modt(uint64_t& x, double a, double scale, uint64_t coeff_modulus);

void print_string_hex(string res);

template <class T>
string to_string(T t, ios_base & (*f)(ios_base&))
{
    ostringstream oss;
    oss << f << t;
    return oss.str();
}

bool compare_strings(vector<string> st1, vector<string> st2);

/*------------------------
 Print the estimated results
 ------------------------*/
void print_labels(string filename, dmat ypred, vector<string> target_geno_id);
void print_estimates(string filename, dmat ypred, vector<string> target_geno_id);
void print_data(dmat yreal, dmat ypred, const string filename = "none");

/*------------------------
 Accuracy measure
 ------------------------*/
void get_dummies(dmat& ylabel, dvec ytest);
void get_labels(dmat& ypred_label, dvec ypred);
void get_mae(double& mae, dvec yreal, dvec ypred);
void get_mse(double& mse, dvec yreal, dvec ypred);

double test_allvariant_macro_acc(dmat ytest, dmat ypred);
double test_nonrefgenotype_macro_acc(dmat ytest, dmat ypred);
void test_mae(double& mean_errs, dvec& errs, dmat yreal, dmat ypred, const string filename = "none");
void test_mse(double& mean_errs, dvec& errs, dmat yreal, dmat ypred, const string filename = "none");

#endif
