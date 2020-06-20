/*
 * @file       HEmpute_bfv.h, header file
 * @brief      defining functions for performing imputation on the bfv scheme
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


class bfvHEmpute
{
public:
    
    Encryptor& encryptor;
    Evaluator& evaluator;
    Decryptor& decryptor;
    
    bfvHEmpute(Encryptor& encryptor, Evaluator& evaluator, Decryptor& decryptor):  encryptor(encryptor), evaluator(evaluator), decryptor(decryptor) {}
    
    /*------------------------
     Imputation for The original data
     ------------------------*/
    
    void encrypt_data(vector<Ciphertext>& encXData, vector<vector<int>> tag_geno_data,
                      vector<long> tag_model_starting_index, long dim);
    
    
    void HEimpute(vector<Ciphertext>& res, vector<Ciphertext> encXData,
                  dvec model0, dmat model, long n_test, double W0scale, double Wscale, size_t plaintext_modulus, vector<long> tag_index);
    
    void decode_Vector(vector<double>& msg, vector<uint64_t> umsg, size_t plaintext_modulus, double scale);
    
    void decrypt_impute(dmat& res, vector<Ciphertext> encres, size_t plaintext_modulus, double scale);

    /*------------------------
     Imputation for stratified population data
     ------------------------*/
    
    void encrypt_parallel(vector<Ciphertext>& encXData, vector<vector<int>> tag_geno_data,
                         vector<long> tag_model_starting_index, long dim, long nparallel);
    
    void HEimpute(vector<Ciphertext>& res, vector<bool>& encres_index,
                  vector<Ciphertext> encXData, dvec model0, dmat model, long n_test,
                  double W0scale, double Wscale, size_t plaintext_modulus, vector<long> tag_index);
    
    void HEimpute_parallel(vector<Ciphertext>& res, vector<bool>& encres_index,
                           vector<Ciphertext> encXData, dvec model0, dmat model, long n_test, long nparallel,
                           double W0scale, double Wscale, size_t plaintext_modulus, vector<long> tag_index);
    
    void decrypt_impute(dmat& res,
                        vector<Ciphertext> encres, vector<bool> encres_index,
                        dvec model0, size_t plaintext_modulus, double scale);
    
};


