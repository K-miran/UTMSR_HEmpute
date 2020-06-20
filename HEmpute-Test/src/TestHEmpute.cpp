/*
 * @file       TestHEmpute.cpp, cpp file
 * @brief      defining functions for preprocessing data
 * @date       Jan. 17, 2020
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"
#include "time.h"
#include <chrono>
#include <sys/resource.h>   // check the memory
#include <random>
#include <iomanip>          // std::setprecision

#include "param.h"
#include "utils.h"
#include "thread.h"
#include "HEmpute_bfv.h"
#include "HEmpute_ckks.h"
#include "TestHEmpute.h"

#define DEBUG false
#define fast true   // using the fast multiplications between single-precision numbers
using namespace std;
 
/*
 Perform the genotype imputation using the bfv scheme on the ALL data
 
 @param[in] model0, [nsnptarget_model] = [p]
 @param[in] model, [p][d]
 @param[in] tag_geno_data, the tag genotype data which will be encrypted using bfv-style
 @param[in] tag_index[nsnptarget], the index that indicates the starting corresponding points in the tag data for each target variant
 this information will only be used for the evaluation phase
 @param[out] ypred, the predicted results of the whole dataset
*/

void TestHEmpute::bfv_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data, vector<long> tag_index){
    
    chrono::high_resolution_clock::time_point time_start, time_end;
    struct rusage usage;
    
    long p = model.size();          // nsnptarget_model = 80,882
    
    // HE parameters
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
    int logq = 27;
    
    // The plain_modulus does not play much of a role;
    // we choose some reasonable value t > results so that it satisfies (results mod t) = results
    size_t plaintext_modulus = (1 << 10);
    
    double Xscale = 2.0;
    double Wscale =  (double)(1 << 6);       // precision of bits (for the parameters Wdata)
    double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg
    
    parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {logq}));
    parms.set_plain_modulus(plaintext_modulus);
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << setprecision(4) << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
#if (DEBUG)
    cout << "(n, logq, logt) = (" << poly_modulus_degree << "," << logq << "," << (long) log2(plaintext_modulus) << ")" << endl;
#endif
    
    HEParam HEparam;
    getHEParam(HEparam, context, DATAparam.dim, DATAparam.n_test, DATAparam.n_snptag);
    
    cout << "+------------------------------------+" << endl;
    cout << "|         Encryption (data)          |" << endl;
    cout << "+------------------------------------+" << endl;
    
    bfvHEmpute hempute(encryptor, evaluator, decryptor);
    
    time_start = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> encXData;
    hempute.encrypt_data(encXData, tag_geno_data, tag_index, DATAparam.dim);

#if (DEBUG)
    cout << "nctxts = "  << encXData.size() << endl;
#endif
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
#if (DEBUG)
    cout << "Modulus chain index for encryption: q[" << context->get_context_data(encXData[0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption: " << decryptor.invariant_noise_budget(encXData[0]) << " bits" << endl;
#endif
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<Ciphertext> encres (p);
   
    hempute.HEimpute(encres, encXData,
                     model0, model, DATAparam.n_test, W0scale, Wscale, plaintext_modulus, tag_index);
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

#if (DEBUG)
    cout << "Modulus chain index for encrypted_res: q[" << context->get_context_data(encres[0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in encrypted_res: " << decryptor.invariant_noise_budget(encres[0]) << " bits" << endl;
#endif
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    hempute.decrypt_impute(ypred, encres, plaintext_modulus, W0scale);

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
}

/*
 Perform the genotype imputation using the ckks scheme on the ALL data
 
 @param[in] model0, [nsnptarget_model] = [p]
 @param[in] model, [p][d]
 @param[in] tag_geno_data, the tag genotype data which will be encrypted using ckks-style
 @param[in] tag_index[nsnptarget], the index that indicates the starting corresponding points in the tag data for each target variant
 this information will only be used for the evaluation phase
 @param[out] ypred, the predicted results of the whole dataset
*/
void TestHEmpute::ckks_HEmpute(dmat& ypred, dvec model0, dmat model,
                               vector<vector<int>> tag_geno_data, vector<long> tag_index){

    chrono::high_resolution_clock::time_point time_start, time_end;
    struct rusage usage;

    long p = model.size();

    // HE parameters
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = (1 << 10);
    int logq = 27;
    size_t plaintext_modulus = (1 << 20); // just for creating the context (no meaning in the ckks scheme)

    double Xscale = (double)(1 << 16);      // tag genotype will be scaled so that the encryption noise doest not affect precision
    double Wscale = (double)(1 << 6);       // precision of bits (for the parameters Wdata)
    double W0scale = Xscale * Wscale;      // scale of the output ciphertext, res/(W0scale) -> msg

    parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {logq}));
    parms.set_plain_modulus(plaintext_modulus);

    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << setprecision(4) << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

    auto &context_data = *context->key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    uint64_t ciphertext_modulus = coeff_modulus[0].value();  // the actual q (coefficient modulus)

    HEParam HEparam;
    getHEParam(HEparam, context, DATAparam.dim, DATAparam.n_test, DATAparam.n_snptag);

    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;

    ckksHEmpute hempute(encryptor, evaluator, decryptor);
    time_start = chrono::high_resolution_clock::now();

    vector<Ciphertext> encXData;
    hempute.encrypt_data(encXData, tag_geno_data, tag_index, DATAparam.dim, Xscale);
    
#if (DEBUG)
    cout << "nctxts = "  << encXData.size() << endl;
#endif
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    vector<Ciphertext> encres (p);
    hempute.HEimpute(encres, encXData,
                     model0, model, DATAparam.n_test, W0scale, Wscale,
                     ciphertext_modulus, tag_index);


    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;

    time_start = chrono::high_resolution_clock::now();

    hempute.decrypt_impute(ypred, encres, ciphertext_modulus, W0scale);

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
}


/*------------------------
 Population stratified data
 ------------------------*/

/*
 Perform the genotype imputation using the bfv scheme on the stratified population data
 
 @param[in] model0, [nsnptarget_model] = [p]
 @param[in] model, [p][d]
 @param[in] tag_geno_data, the tag genotype data which will be encrypted using bfv-style
 @param[in] tag_index[nsnptarget], the index that indicates the starting corresponding points in the tag data for each target variant
 this information will only be used for the evaluation phase
 @param[out] ypred, the predicted results of the whole dataset
 */
void TestPOPHEmpute::bfv_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data,
                                  vector<long> tag_index, bool parallel)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    struct rusage usage;
    
    long p = model.size();
    
    // HE parameters
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = (1 << 10);  // max coeff_modulus bit-lenght =  109
    int logq = 27;
    
    // In this example the plain_modulus does not play much of a role; we choose some reasonable value
    // res < t in order to satisfy (res mod t) = res
    size_t plaintext_modulus = (1 << 10);
    
    double Xscale = 2.0;
    double Wscale =  (double)(1 << 6);       // precision of bits (for Wdata)
    double W0scale = Xscale * Wscale;       // scale of the output ciphertext, res/(W0scale) -> msg
    
    parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {logq}));
    parms.set_plain_modulus(plaintext_modulus);
    
    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << setprecision(4) << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
#if (DEBUG)
    cout << "(n, logq, logt) = (" << poly_modulus_degree << "," << logq << "," << (long) log2(plaintext_modulus) << ")" << endl;
#endif
    
    HEParam HEparam;
    getHEParam(HEparam, context, DATAparam.dim, DATAparam.n_test, DATAparam.n_snptag);
    long nparallel = floor((double)(poly_modulus_degree)/(double)(DATAparam.n_test));
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    bfvHEmpute hempute(encryptor, evaluator, decryptor);
    time_start = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> encXData;
    if (parallel){
        hempute.encrypt_parallel(encXData, tag_geno_data, tag_index, DATAparam.dim, nparallel);
    }
    else{
        hempute.encrypt_data(encXData, tag_geno_data, tag_index, DATAparam.dim);
    }
    
#if (DEBUG)
    cout << "nctxts = "  << encXData.size() << endl;
#endif
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
#if (DEBUG)
    cout << "Modulus chain index for encryption: q[" << context->get_context_data(encXData[0].parms_id()) ->chain_index() << "]" << endl;
    cout << "Noise budget in fresh encryption: " << decryptor.invariant_noise_budget(encXData[0]) << " bits" << endl;
#endif
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> encres (p);
    vector<bool> encres_index (p, false);
    
    if(parallel){
        hempute.HEimpute_parallel(encres, encres_index,
                                  encXData, model0, model,
                                  DATAparam.n_test, nparallel, W0scale, Wscale,
                                  plaintext_modulus, tag_index);
    }
    else{
        hempute.HEimpute(encres, encres_index,
                         encXData, model0, model,
                         DATAparam.n_test, W0scale, Wscale,
                         plaintext_modulus, tag_index);
    }
    
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    hempute.decrypt_impute(ypred, encres, encres_index,
                           model0, plaintext_modulus, W0scale);
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
}


/*
 Perform the genotype imputation using the ckks scheme on the stratified population data
 
 @param[in] model0, [nsnptarget_model] = [p]
 @param[in] model, [p][d]
 @param[in] tag_geno_data, the tag genotype data which will be encrypted using ckks-style
 @param[in] tag_index[nsnptarget], the index that indicates the starting corresponding points in the tag data for each target variant
 this information will only be used for the evaluation phase
 @param[out] ypred, the predicted results of the whole dataset
 */
void TestPOPHEmpute::ckks_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data,
                                  vector<long> tag_index, bool parallel)
{
    chrono::high_resolution_clock::time_point time_start, time_end;
    struct rusage usage;
    
    long p = model.size();
    
    // HE parameters
    EncryptionParameters parms(scheme_type::BFV);
    size_t poly_modulus_degree = (1 << 10);
    int logq = 27;
    size_t plaintext_modulus = (1 << 20); // just for creating the context (no meaning in the ckks scheme)
    
    double Xscale = (double)(1 << 16);      // tag genotype will be scaled so that the encryption noise doest not affect precision
    double Wscale = (double)(1 << 6);       // precision of bits (for Wdata)
    double W0scale = Xscale * Wscale;      // scale of the output ciphertext, res/(W0scale) -> msg
    
    parms.set_poly_modulus_degree(poly_modulus_degree); // n = degree
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {logq}));
    parms.set_plain_modulus(plaintext_modulus);

    cout << "+------------------------------------+" << endl;
    cout << "|           Key Generation           |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    auto context = SEALContext::Create(parms);
    KeyGenerator keygen(context);
    auto public_key = keygen.public_key();
    auto secret_key = keygen.secret_key();
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Scheme generation (milliseconds) : " << time_diff.count()/1000.0 << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << setprecision(4) << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
    auto &context_data = *context->key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    uint64_t ciphertext_modulus = coeff_modulus[0].value();  // the actual q (coefficient modulus)
    
    HEParam HEparam;
    getHEParam(HEparam, context, DATAparam.dim, DATAparam.n_test, DATAparam.n_snptag);
    long nparallel = floor((double)(poly_modulus_degree)/(double)(DATAparam.n_test));
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Encryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    ckksHEmpute hempute(encryptor, evaluator, decryptor);
    time_start = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> encXData;
    if (parallel){
        hempute.encrypt_parallel(encXData, tag_geno_data, tag_index, DATAparam.dim, nparallel, Xscale);
    }
    else{
        hempute.encrypt_data(encXData, tag_geno_data, tag_index, DATAparam.dim, Xscale);
    }
    
#if (DEBUG)
    cout << "nctxts = "  << encXData.size() << endl;
#endif
    
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Encryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
    
    cout << "+------------------------------------+" << endl;
    cout << "|             Evaluation             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    vector<Ciphertext> encres (p);
    vector<bool> encres_index (p, false);
    
    if(parallel){
        hempute.HEimpute_parallel(encres, encres_index,
                                  encXData, model0, model,
                                  DATAparam.n_test, nparallel, W0scale, Wscale,
                                  ciphertext_modulus, tag_index);
    }
    else{
        hempute.HEimpute(encres, encres_index,
                         encXData, model0, model,
                         DATAparam.n_test, W0scale, Wscale,
                         ciphertext_modulus, tag_index);
    }
 
 
    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Evaluation (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
 
    cout << "+------------------------------------+" << endl;
    cout << "|             Decryption             |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    hempute.decrypt_impute(ypred, encres, encres_index,
                           model0, ciphertext_modulus, W0scale);

    time_end = chrono::high_resolution_clock::now();
    time_diff = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "Decryption (seconds) : " << time_diff.count()/(1000000.0) << endl;
    getrusage(RUSAGE_SELF, &usage);
    cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
}
