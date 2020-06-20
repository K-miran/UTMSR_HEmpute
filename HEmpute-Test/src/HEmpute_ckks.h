/*
 * @file       HEmpute_ckks.h, header file
 * @brief      defining functions for performing imputation on the ckks scheme
 *
 * @author     Miran Kim
 * @date       June. 4, 2019
 * @copyright  GNU Pub License
 */

#include <vector>
#include <complex>
#include <iostream>
#include <string>
#include <cstddef>

#include "utils.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;


class ckksHEmpute
{
public:
    
    Encryptor& encryptor;
    Evaluator& evaluator;
    Decryptor& decryptor;
    
    ckksHEmpute(Encryptor& encryptor, Evaluator& evaluator, Decryptor& decryptor): encryptor(encryptor), evaluator(evaluator), decryptor(decryptor) {}
    
    /*------------------------
     Imputation for The original data
     ------------------------*/
    
    void encrypt_data(vector<Ciphertext>& encXData, vector<vector<int>> tag_geno_data,
                      vector<long> tag_model_starting_index, long dim, double Xscale);
    
    void HEimpute(vector<Ciphertext>& res, vector<Ciphertext> encXData,
                  dvec model0, dmat model, long n_test, double W0scale, double Wscale,
                  uint64_t coeff_modulus, vector<long> tag_index);
    
    void decode_Vector(vector<double>& msg, vector<uint64_t> umsg, uint64_t coeff_modulus, double scale);
    
    void decrypt_impute(dmat& res, vector<Ciphertext> encres, uint64_t coeff_modulus, double scale);

    /*------------------------
     Imputation for stratified population data
     ------------------------*/
    
    void encrypt_parallel(vector<Ciphertext>& encXData, vector<vector<int>> tag_geno_data,
                         vector<long> tag_model_starting_index, long dim, long nparallel, double Xscale);
    
    void HEimpute(vector<Ciphertext>& res, vector<bool>& encres_index,
                  vector<Ciphertext> encXData, dvec model0, dmat model, long n_test, 
                  double W0scale, double Wscale, uint64_t coeff_modulus, vector<long> tag_index);
    
    void HEimpute_parallel(vector<Ciphertext>& res, vector<bool>& encres_index,
                      vector<Ciphertext> encXData, dvec model0, dmat model, long n_test, long nparallel,
                      double W0scale, double Wscale, uint64_t coeff_modulus, vector<long> tag_index);
    
    void decrypt_impute(dmat& res, vector<Ciphertext> encres, vector<bool> encres_index,
                        dvec model0, size_t coeff_modulus, double scale);
    
};
