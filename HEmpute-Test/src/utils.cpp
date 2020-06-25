/*
 * @file       utils.cpp, cpp file
 * @brief      defining functions for reading and storing data
 *
 * @author     Miran Kim
 * @date       June.3, 2019
 * @copyright  GNU Pub License
 */

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <vector>

#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>

#include "utils.h"
#include "thread.h"

using namespace std;

/*
 Return a positive integer x such that
 x <- round(a * scale) mod modulus
 
 @param[in] a, the double-type input value
 @param[in] scale, the scale factor
 @param[in] modulus, the modulus (plaintext modulus in the BFV scheme, or ciphertext modulus in the CKKS scheme)
 @param[out] x, the uint64_t value such that
 if a >= 0, then x = (a * scale)
 if a < 0, then x = (a * scale) coeff_modulus = coeff_mod - (a * scale)
*/
bool scale_and_modt(uint64_t& x, double a, double scale, uint64_t modulus)
{
    if(a == 0){
        x = 0;
        return false;
    }
    else if(a > 0){
        x = (uint64_t) round(a * scale);
        if (x != 0){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        x = (uint64_t) round(- a * scale);
        if(x != 0){
            x = (uint64_t) modulus - x;
            return true;
        }
        else{
            return false;
        }
    }
    
//    uint64_t tmp = (uint64_t) round(a * scale);
//    x = tmp - (coeff_modulus & static_cast<uint64_t>(-static_cast<int64_t>(tmp >= coeff_modulus)));
}


void print_string_hex(string res)
{
    cout << hex;
    int size = res.size();
    for (int i= 0; i < size; i++){
        cout << (int) res[i] << " ";
    }
    cout << endl;
}

bool compare_strings(vector<string> st1, vector<string> st2){
    bool res;
    size_t len = st1.size();
    size_t i = 0;
    
    while(1){
        if(st1[i] != st2[i]){
            res = false;
            break;
        }
        i++;
        if(i == len){
            res = true;
            break;
        }
    }
    return res;
}


/*
 Print the predicted labels from the predicted results.
 
 @param[in] ypred, the predicted real-valued array, [p][n]
 @param[in] target_geno_id, the name of the target SNPs
 @param[out] ypred_label
 */
void print_labels(string filename, dmat ypred, vector<string> target_geno_id)
{
    long nsnp = ypred.size();
    long ntest = ypred[0].size();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file
    
    outf << "Subject ID,target SNP,0,1,2" << endl;
    for(long i = 0; i < ntest; ++i){
        vector<vector<bool>> ypred_label(nsnp, vector<bool> (3));
        
        MT_EXEC_RANGE(nsnp, first, last);
        for(long j = first; j < last; ++j){
            double pred_value = ypred[j][i];
            if(pred_value <= 0.25){
                ypred_label[j] = {true, false, false};
            }
            else if(pred_value < 0.75){
                ypred_label[j] = {false, true, false};
            }
            else {
                ypred_label[j] = {false, false, true};
            }
        }
        MT_EXEC_RANGE_END;
        
        for(long j = 0; j < nsnp; ++j){
            outf << i << "," << target_geno_id[j] << ",";
            outf << ypred_label[j][0] << "," << ypred_label[j][1] << "," << ypred_label[j][2] << endl;
        }
    }

    outf.close();
}

/*
 @param[in] ypred, the predicted real-valued array, [p][n]
 @param[out] ypred_label s.t.
 print the estimated values 
 */
void print_estimates(string filename, dmat ypred, vector<string> target_geno_id)
{
    long nsnp = ypred.size();
    long ntest = ypred[0].size();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file
    
    outf << "Subject ID,target SNP,estimates" << endl;
    
    for(long i = 0; i < ntest; ++i){
        for(long j = 0; j < nsnp; ++j){
            outf << i << "," << target_geno_id[j] << "," << ypred[j][i] << endl;
        }
    }
    
    outf.close();
}

/*
 Print the actual genotypes and the predicted values.
 
 @param[in] yreal, the actual real-valued array, [p][n]
 @param[in] ypred, the predicted real-valued array, [p][n]
 @param[in] filename, the filename
 */
void print_data(dmat yreal, dmat ypred, string filename){
    long nsnp = yreal.size();
    
    fstream outf;
    outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);
    for(long i = 0; i < nsnp; ++i){
        long ntest = yreal[i].size();
        for(long j = 0; j < ntest; ++j){
            outf << yreal[i][j] << "\t" << ypred[i][j] << endl;
        }
    }
}

/*
 Change the actual vector into class vectors
 
 @param[in] ytest, the double-type actual vector
 @param[out] ylabel[n][3], the data label in a way that
 ylabel[i] = (1,0,0) when ytest[i] = 0;
 ylabel[i] = (0,1,0) when ytest[i] = 0.5;
 ylabel[i] = (0,0,1) when ytest[i] = 1.0.
 */
void get_dummies(dmat& ylabel, dvec ytest){
    long ntest = ytest.size();
    
    for(long i = 0; i < ntest; ++i){
        dvec tmp;
        if(ytest[i] == 0.0){
            tmp.push_back(1.0);
            tmp.push_back(0.0);
            tmp.push_back(0.0);
        }
        else if(ytest[i] == 0.5){
            tmp.push_back(0.0);
            tmp.push_back(1.0);
            tmp.push_back(0.0);
        }
        else {
            tmp.push_back(0.0);
            tmp.push_back(0.0);
            tmp.push_back(1.0);
        }
        ylabel.push_back(tmp);
    }
}

/*
 @param[in] ypred, the predicted real-valued vector
 @param[out] ypred_label s.t. ypred_label[i][j] = sigmoid(1 - |ypred[i]-(j/2)|)
 */
void get_labels(dmat& ypred_label, dvec ypred)
{
    long npred = ypred.size();
    double val_list[3] = {0.0, 0.5, 1.0};
    
    for(long i = 0; i < npred; ++i){
        dvec dtmp;
        for(long j = 0; j < 3; ++j){
            double diff = abs(ypred[i] - val_list[j]);
            double dval = 1.0 - diff;
            dtmp.push_back(1.0 / (1.0 + exp(-dval)));
        }
        ypred_label.push_back(dtmp);
    }
}


/*
 Calculate the mean absolute errors of the current SNP variant
 
 @param[in] yreal, the actual genotypes
 @param[in] ypred, the predicted genotypes
 @param[out] mae, mean absolute errors of two input vectors (excludiing the cases that yreal is zero)
 @throws std::invalid_argument if the sizes of yreal and yprea do not match
 */
void get_mae(double& mae, dvec yreal, dvec ypred){
    long ntest = yreal.size();
    if(ntest != (long) ypred.size()){
        throw std::invalid_argument("mae: dimension mismatch");
    }
    else{
        long num_non_ref_geno_indices = 0;
        mae = 0.0;
        for(long i = 0; i < ntest; ++i){
            //if(yreal[i] > 0.0){
                double rtmp = abs(yreal[i] - ypred[i]);
                //cout << i << ": " << ytest[i] << "," << ypred[i] << endl;
                mae += rtmp;
                num_non_ref_geno_indices++;
            //}
        }
        mae /= (double)(num_non_ref_geno_indices);
    }
}

/*
 Calculate the mean absolute errors of the current SNP variant
 
 @param[in] yreal, the actual genotypes
 @param[in] ypred, the predicted genotypes
 @param[out] mse, mean squared  errors of two input vectors (excludiing the cases that yreal is zero)
 @throws std::invalid_argument if the sizes of yreal and yprea do not match
 */
void get_mse(double& mse, dvec yreal, dvec ypred){
    long ntest = yreal.size();
    if(ntest != (long) ypred.size()){
        throw std::invalid_argument("mse: dimension mismatch");
    }
    else{
        long num_non_ref_geno_indices = 0;
        mse = 0.0;
        for(long i = 0; i < ntest; ++i){
            if(yreal[i] > 0.0){
                double rtmp = pow(yreal[i] - ypred[i], 2.0);
                mse += rtmp;
                num_non_ref_geno_indices++;
            }
        }
        mse /= (double)(num_non_ref_geno_indices);
    }
}


/*
 Calculate the macro-accuracy of all the SNP variants
 by computing the average of acc[j]'s
 s.t. acc[j] = [#(correctly imputed individuals for variant j)/#(individuals for variant j)]
 
 @param[in] yreal, the actual genotypes of size [p][n]
 @param[in] ypred, the predicted genotypes of size [p][n]
 @param[out] res, the macro-aggregated accuracy
 */
double test_allvariant_macro_acc(dmat yreal, dmat ypred)
{
    long nsnp = yreal.size();
    vector<double> var_acc (nsnp);
    
    // 1. compute each variant accuracy
    // #(correctly imputed individuals for variant j)/#(individuals for variant j)
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        long n_match = 0;
        long nsamples = yreal[j].size();
        
        for(long i = 0; i < nsamples; ++i){
            if(yreal[j][i] == 0.0){
                if((ypred[j][i] < 0.25)){
                    n_match++;
                }
            }
            else if(yreal[j][i] == 0.5){
                if((ypred[j][i] >= 0.25) && (ypred[j][i] <= 0.75)){
                    n_match++;
                }
            }
            else if(yreal[j][i] == 1.0){
                if(ypred[j][i] >= 0.75){
                    n_match++;
                }
            }
        }
        var_acc[j] = (double)(n_match)/ (double)(nsamples);
    }
    MT_EXEC_RANGE_END;
    
    // 2. Take the average of the variant accuracies
    double res = 0.0;
    for(long j = 0; j < nsnp; ++j){
        res += var_acc[j];
    }
    res /= (double) nsnp;
    return res;
}

/*
 Calculate the macro-accuracy of non-reference genotypes
 by computing the average of acc[j]'s
 s.t. acc[j] = [#(correctly imputed non.ref individuals for variant j)/#(non.ref individuals for variant j)]
 
 @param[in] yreal, the actual genotypes of size [p][n]
 @param[in] ypred, the predicted genotypes of size [p][n]
 @param[out] res, the macro-aggregated accuracy
 */
double test_nonrefgenotype_macro_acc(dmat yreal, dmat ypred)
{
    long nsnp = yreal.size();
    vector<long> n_non_ref (nsnp, 0ULL);  // #(non.ref individuals for variant j)
    vector<long> n_non_ref_match (nsnp, 0ULL);  // #(correctly imputed non.ref individuals for variant j)
    vector<double> var_acc (nsnp);
    
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        long ntest = yreal[j].size();
        for(long i = 0; i < ntest; ++i){
            if(yreal[j][i] == 0.5){
                n_non_ref[j]++;
                if((ypred[j][i] >= 0.25) && (ypred[j][i] <= 0.75)){
                    n_non_ref_match[j]++;
                }
            }
            else if(yreal[j][i] == 1.0){
                n_non_ref[j]++;
                if(ypred[j][i] >= 0.75){
                    n_non_ref_match[j]++;
                }
            }
        }
        if(n_non_ref[j] != 0){
            var_acc[j] = (double)(n_non_ref_match[j])/ (double)(n_non_ref[j]);
        }
    }
    MT_EXEC_RANGE_END;
    
#if 0
    for(long j = 0; j < 10; ++j){
        cout << var_acc[j] << "=" << n_non_ref_match[j] << "/" << n_non_ref[j] << endl;
    }
#endif
    
    // 2. Take the average over the variant accuracies
    double res = 0.0;
    long n_pred = 0;
    for(long j = 0; j < nsnp; ++j){
        if(n_non_ref[j] != 0){
            res += var_acc[j];
            n_pred++;
        }
    }
    res /= (double) n_pred;
    return res;
}

/*
 Calculate the mean absolute errors of all the SNP variants
 as well as their average.
 
 @param[in] yreal, the actual genotypes of size [p][n]
 @param[in] ypred, the predicted genotypes of size [p][n]
 @param[in] filename, the filename where mAUCs are stored
 @params[out] mean_errs, the average of the mean absolute errors of all the variants
 @param[out] errs, the mean absolute errors of all the variants
 */
void test_mae(double& mean_errs, dvec& errs, dmat yreal, dmat ypred, string filename){
    long nsnp = yreal.size();
    mean_errs = 0.0;

    MT_EXEC_RANGE(nsnp, first, last);
    for(long i = first; i < last; ++i){
        get_mae(errs[i], yreal[i], ypred[i]);
    }
    MT_EXEC_RANGE_END;
    
    if(filename == "none"){
        for(long i = 0; i < nsnp; ++i){
            mean_errs = ((double) i/(i+1)) * mean_errs + (1.0/(i+1)) * errs[i];
        }
    }
    else{
        fstream outf;
        outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file

        for(long i = 0; i < nsnp; ++i){
            outf << i << ": " << errs[i] << endl;
            mean_errs = ((double) i/(i+1)) * mean_errs + (1.0/(i+1)) * errs[i];
        }
        outf.close();
    }
}

/*
 Calculate the mean squared errors of all the SNP variants
 as well as their average.
 
 @param[in] yreal, the actual genotypes of size [p][n]
 @param[in] ypred, the predicted genotypes of size [p][n]
 @param[in] filename, the filename where the results are stored
 @params[out] mean_errs, the average of the mean squared errors of all the variants
 @param[out] errs, the mean squared errors of all the variants
 */
void test_mse(double& mean_errs, dvec& errs, dmat yreal, dmat ypred, string filename){
    long nsnp = yreal.size();
    mean_errs = 0.0;

    MT_EXEC_RANGE(nsnp, first, last);
    for(long i = first; i < last; ++i){
        get_mse(errs[i], yreal[i], ypred[i]);
    }
    MT_EXEC_RANGE_END;
    
    if(filename == "none"){
        for(long i = 0; i < nsnp; ++i){
            mean_errs = ((double) i/(i+1)) * mean_errs + (1.0/(i+1)) * errs[i];
        }
    }
    else{
        fstream outf;
        outf.open(filename.c_str(), fstream::in | fstream::out | fstream::app);   // open the file

        for(long i = 0; i < nsnp; ++i){
            outf << i << ": " << errs[i] << endl;
            mean_errs = ((double) i/(i+1)) * mean_errs + (1.0/(i+1)) * errs[i];
        }
        outf.close();
    }
}
