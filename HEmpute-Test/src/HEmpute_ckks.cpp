/*
 * @file       HEmpute.cpp, cpp file
 * @brief      defining functions for performing imputation on the ckks scheme
 * @date       June. 4, 2019
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
#include <sys/resource.h>
#include <random>
#include <iomanip>

#include "utils.h"
#include "thread.h"
#include "HEmpute_ckks.h"

#define DEBUG false
using namespace std;

/*
 Encrypt tagSNPs using tag_model_starting_index information.
 
 @param[in] tag_geno_data, [p][n] = [16184][1004] and each entry is an elemenet of {0, 1, 2}
 @param[in] tag_model_starting_index, starting index of tagSNPs
 @param[in] dim, the number of predictors
 @params[in] Xscale, the scaling factor of tagSNPs
 @param[out] encXData[p'], the encryptions of the tag dataset
 in a way that each column of the snp variant is encrypted as a single ciphertext
 */
void ckksHEmpute::encrypt_data(vector<Ciphertext>& encXData,
                               vector<vector<int>> tag_geno_data, vector<long> tag_model_starting_index, long dim, double Xscale)
{
    long nshift = static_cast<long> (log2(Xscale) - 1.0);
    long length = tag_model_starting_index.size();
    long nsnp = (tag_model_starting_index[length - 1] + dim);   // the last starting one plus dimension (the actual number of tag genotypes for imputation)

    // Each entry of tag_geno_data is an element of {0,1,2} and we want to scale it by a factor of Xscale
    // But the parameters were trained on the data with the values of {0,0.5,1},
    // so we scale the data by the half.
    encXData.resize(nsnp);
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        encryptor.encrypt_ckks(tag_geno_data[j], nshift, encXData[j], MemoryPoolHandle().New(false));
    }
    MT_EXEC_RANGE_END;
}

/*
 Perform the genotype imputation over encrypted data
 It requires at most (dim) constant multiplications per SNP.
 This implies that the total complexity is (n_snptarget) * (dim)
 Note that it takes around 16 microsecond to perform one scalar multiplicaiton
 
 @param[in] encXData[ptag], the encryption of tagSNPs
 @param[in] model0[p], the intercepts where p denotes the number of the actual imputed genotypes (p <= 80882)
 @param[in] model[p][d], the trainded parameters where d denotes the dimension of the actual trained models (d = 2 * vicinities)
 @param[in] n_test, the number of test samples
 @param[in] W0scale, the scaling factor of the intercept w0
 @param[in] Wscale, the scaling factor of the model parameterss
 @param[in] coeff_modulus, the ciphertext modulus (q) in the ckks scheme
 @param[in] tag_index[p], the index that indicates the starting corresponding points in the tag data for each target variant
 @param[out] res[p], the predicted encrypted results for each snp variant.
 */
void ckksHEmpute::HEimpute(vector<Ciphertext>& res, vector<Ciphertext> encXData,
                           dvec model0, dmat model, long n_test, double W0scale, double Wscale,
                           uint64_t coeff_modulus,  vector<long> tag_index)
{
    long nsnp = model.size();       // the number of the actual imputed genotypes
    long dim = model[0].size();     // the dimension of the actual trained models

   
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        // res[j] <- model[j][0] * encXdata[j_tag] + ... + model[j][dim-1] * encXdata[j_tag+dim-1],
        long k_start = 0;
        
        // find the first nonzero weight
        uint64_t coeff;
        while(1){
            if(scale_and_modt(coeff, model[j][k_start], Wscale, coeff_modulus)){
                break;
            }
            k_start++;
        }
        long j_tag = tag_index[j];
        res[j] = encXData[j_tag + k_start];
        evaluator.multiply_scalar_inplace_ckks_fast(res[j], coeff, MemoryPoolHandle().New(false));
        
        // add the remaining terms
        long k = k_start + 1;
        while(k < dim){
            if(scale_and_modt(coeff, model[j][k], Wscale, coeff_modulus)){
                Ciphertext ctemp = encXData[j_tag + k];
                evaluator.multiply_scalar_inplace_ckks_fast(ctemp, coeff, MemoryPoolHandle().New(false));
                evaluator.add_inplace(res[j], ctemp);
            }
            k++;
        }
        
        // add the (nonzero) intercepts by computing
        // res[j] += sum_{0<=i<ntest} coeff * x^i
        if(scale_and_modt(coeff, model0[j], W0scale, coeff_modulus)){
            evaluator.add_plain_inplace_ckks(res[j], coeff, n_test);
        }
    }
    MT_EXEC_RANGE_END;
}

/*
 Divide the input vector by plaintext_modulus and apply the modular reduction.
 
 @param[in] umsg, the vector (obtained after the decryption of the bfv scheme)
 @param[in] coeff_modulus, the coefficient modulus in the ckks scheme
 @param[in] scale, the scaling factor of the messages
 @param[out] msg, double-type vector such that msg[i] = umsg[i]/(scale) mod plaintext_modulus
 */
void ckksHEmpute::decode_Vector(vector<double>& msg, vector<uint64_t> umsg, uint64_t coeff_modulus, double scale){
    long n = umsg.size();
    double half_coeff_modulus = (coeff_modulus + 1) >> 1;
    
    for(long i = 0; i < n; ++i){
        if(umsg[i] < half_coeff_modulus){
            double tmp = static_cast<double> (umsg[i]) / scale;
            msg.push_back(tmp);
        }
        else{
            uint64_t utmp = coeff_modulus - umsg[i];
            double tmp = static_cast<double> (utmp) / scale;
            msg.push_back(-tmp);
        }
    }
}

/*
 Decrypt the input ciphertexts and then decode the resulting vectors.
 
 @param[in] encres[p], the encrypted results
 @param[in] plaintext_modulus, the plaintext modulus
 @param[in] scale, the scaling factor of the messages
 @param[out] res, the decrypted results
 */
void ckksHEmpute::decrypt_impute(dmat& res, vector<Ciphertext> encres, uint64_t coeff_modulus, double scale)
{
    long nsnp = res.size();
    long n_test = res[0].size();
    
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        Plaintext plaintmp;
        decryptor.decrypt_ckks(encres[j], plaintmp);
        uint64_t* plainvec = plaintmp.data();       // convert hex into decimal
        
        vector<uint64_t> umsg;
        for(long i = 0; i < n_test; ++i){
            umsg.push_back(plainvec[i]);
        }
        vector<double> dmsg;
        decode_Vector(dmsg, umsg, coeff_modulus, scale);
        res[j] = dmsg;
    }
    MT_EXEC_RANGE_END;
}


/*
 Encrypt tagSNPs using tag_model_starting_index information.
 
 @param[in] tag_geno_data, [p][n] = [16184][1004] and each entry is an elemenet of {0, 1, 2}
 @param[in] tag_model_starting_index, starting index of tagSNPs
 @param[in] dim, the number of predictors
 @param[in] nparallel, the number of SNP variants to be encrypted into a single ciphertext
 @param[in] Xscale, the scaling factor of the input
 @param[out] encXData[p'], the encryptions of the tag dataset
 */
void ckksHEmpute::encrypt_parallel(vector<Ciphertext>& encXData, vector<vector<int>> tag_geno_data,
                                  vector<long> tag_model_starting_index, long dim, long nparallel, double Xscale)
{
    long nshift = static_cast<long> (log2(Xscale) - 1.0);
    long length = tag_model_starting_index.size();
    long nsnp = (tag_model_starting_index[length - 1] + dim);   // the last starting one plus dimension (the actual number of tag genotypes for imputation)
   
    long nctxts = ceil((double)(nsnp)/(double)(nparallel));
    long nparallel_last = nsnp % nparallel;
    if(nparallel_last == 0){
        nparallel_last = nparallel;
    }
    
    // divide the whole dataset and encrypt several columns together (as many as possible)
    encXData.resize(nctxts);
    
    MT_EXEC_RANGE(nctxts - 1, first, last);
    for(long j = first; j < last; ++j){
        vector<int> merged_tag_geno_data;
        long jst = j * nparallel;
        for(long i = 0; i < nparallel; ++i){
            merged_tag_geno_data.insert(merged_tag_geno_data.end(),
                                        tag_geno_data[i + jst].begin(), tag_geno_data[i + jst].end());
        }
        encryptor.encrypt_ckks(merged_tag_geno_data, nshift,
                               encXData[j], MemoryPoolHandle().New(false));
    }
    MT_EXEC_RANGE_END;
    
    vector<int> merged_tag_geno_data;
    long jst = (nctxts - 1) * nparallel;
    for(long i = 0; i < nparallel_last; ++i){
        merged_tag_geno_data.insert(merged_tag_geno_data.end(),
                                    tag_geno_data[i + jst].begin(), tag_geno_data[i + jst].end());
    }
    encryptor.encrypt_ckks(merged_tag_geno_data, nshift,
                           encXData[nctxts - 1], MemoryPoolHandle().New(false));
}

/*
 Perform the genotype imputation over encrypted data
 It requires at most (dim) constant multiplications per SNP.
 This implies that the total complexity is (n_snptarget) * (dim)
 Note that it takes around 16 microsecond to perform one scalar multiplicaiton
 
 @param[in] encXData[ptag], the encryption of tagSNPs
 @param[in] model0[p], the intercepts where p denotes the number of the actual imputed genotypes (p <= 80882)
 @param[in] model[p][d], the trainded parameters where d denotes the dimension of the actual trained models (d = 2 * vicinities)
 @param[in] n_test, the number of test sampels
 @param[in] W0scale, the scaling factor of the intercept w0
 @param[in] Wscale, the scaling factor of the model parameterss
 @param[in] coeff_modulus, the coefficient modulus (q) in the ckks scheme
 @param[in] tag_index[p], the index that indicates the starting corresponding points in the tag data for each target variant
 @param[out] res[p], the predicted encrypted results for each snp variant.
 @param[out] encres_index[p], the indicator that represents whether the corresponding tag ciphertexts are used or not for performing the imputation at each target variant.
 Each entry is 1 if at least one scalar multiplication is performed between tagSNP and params; otherwise 0.
 */
void ckksHEmpute::HEimpute(vector<Ciphertext>& res, vector<bool>& encres_index, 
                           vector<Ciphertext> encXData, dvec model0, dmat model, long n_test,
                           double W0scale, double Wscale, uint64_t coeff_modulus, vector<long> tag_index)
{
    long nsnp = model.size();   // the number of the actual imputed genotypes
    long dim = model[0].size(); // the dimension of the actual trained models
    
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        // res[j] <- model[j][0] * encXdata[j_tag] + ... + model[j][dim-1] * encXdata[j_tag+dim-1],
        long k_start = 0;
        uint64_t coeff;
        while(1){
            if(scale_and_modt(coeff, model[j][k_start], Wscale, coeff_modulus)){    // find the first nonzero weight
                break;
            }
            k_start++;
            if(k_start == dim)
                break;
        }
        
        if(k_start != dim) {
            long j_tag = tag_index[j];
            res[j] = encXData[j_tag + k_start];
            evaluator.multiply_scalar_inplace_ckks_fast(res[j], coeff, MemoryPoolHandle().New(false));
            encres_index[j] = true;
            
            // add the remaining terms
            long k = k_start + 1;
            while(k < dim){
                if(scale_and_modt(coeff, model[j][k], Wscale, coeff_modulus)){
                    Ciphertext ctemp = encXData[j_tag + k];
                    evaluator.multiply_scalar_inplace_ckks_fast(ctemp, coeff, MemoryPoolHandle().New(false));
                    evaluator.add_inplace(res[j], ctemp);
                }
                k++;
            }
            
            // add the (nonzero) intercepts by computing
            // res[j] += sum_{0<=i<ntest} coeff * x^i
            if(scale_and_modt(coeff, model0[j], W0scale, coeff_modulus)){
                evaluator.add_plain_inplace_ckks(res[j], coeff, n_test);
            }
        }
    }
    MT_EXEC_RANGE_END;
}

/*
 Perform the genotype imputation over encrypted data
 It requires at most (dim) constant multiplications per SNP.
 This implies that the total complexity is (n_snptarget) * (dim)
 Note that it takes around 16 microsecond to perform one scalar multiplicaiton
 
 @param[in] encXData[ptag], the encryption of tagSNPs
 @param[in] model0[p], the intercepts where p denotes the number of the actual imputed genotypes (p <= 80882)
 @param[in] model[p][d], the trainded parameters where d denotes the dimension of the actual trained models (d = 2 * vicinities)
 @param[in] n_test, the number of test sampels
 @param[in] W0scale, the scaling factor of the intercept w0
 @param[in] Wscale, the scaling factor of the model parameterss
 @param[in] plaintext_modulus, the plaintext modulus (t) in the bfv scheme
 @param[in] tag_index[p], the index that indicates the starting corresponding points in the tag data for each target variant
 @param[out] res[p], the predicted encrypted results for each snp variant.
 @param[out] encres_index[p], the indicator that represents whether the corresponding tag ciphertexts are used or not for performing the imputation at each target variant.
 Each entry is 1 if at least one scalar multiplication is performed between tagSNP and params; otherwise 0.
 */
void ckksHEmpute::HEimpute_parallel(vector<Ciphertext>& res, vector<bool>& encres_index,
                                vector<Ciphertext> encXData, dvec model0, dmat model, long n_test, long nparallel,
                                double W0scale, double Wscale, uint64_t coeff_modulus, vector<long> tag_index)
{
    long nsnp = model.size();   // the number of the actual imputed genotypes
    long dim = model[0].size(); // the dimension of the actual trained models
    cout << "> (n_test, nparallel) = (" << n_test << ", " << nparallel << ")" << endl;
    
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        // res[j] <- model[j][0] * encXdata[j_tag] + ... + model[j][dim-1] * encXdata[j_tag+dim-1],
        long k_start = 0;
        uint64_t coeff;
        while(1){
            if(scale_and_modt(coeff, model[j][k_start], Wscale, coeff_modulus)){    // find the first nonzero weight
                break;
            }
            k_start++;
            if(k_start == dim)
                break;
        }
       
        if(k_start != dim){
            vector<Ciphertext> res_chunk (nparallel);
            vector<bool> indicator (nparallel, false);
            
            long index = tag_index[j] + k_start;
            long i = (long) floor((double)(index)/(double)(nparallel));
            long l = index % nparallel;  // this index shows which entry in the ciphertext we use
            
            res_chunk[l] = encXData[i];
            evaluator.multiply_scalar_inplace_ckks_fast(res_chunk[l], coeff, MemoryPoolHandle().New(false));
            encres_index[j] = true;
            indicator[l] = true;
            
            // add the remaining terms
            long k = k_start + 1;
            while(k < dim){
                if(scale_and_modt(coeff, model[j][k], Wscale, coeff_modulus)){
                    long index = tag_index[j] + k;
                    long i = (long) floor((double)(index)/(double)(nparallel));
                    long l = index % nparallel;
                    
                    // first perform the scalar multiplication with the model
                    // if l-th ciphertext is occupied, then add the output ciphertext to it; otherwise, copy it.
                    Ciphertext ctemp = encXData[i];
                    evaluator.multiply_scalar_inplace_ckks_fast(ctemp, coeff, MemoryPoolHandle().New(false));
                    if(indicator[l]){
                        evaluator.add_inplace(res_chunk[l], ctemp);
                    }
                    else{
                        res_chunk[l] = ctemp;
                    }
                    indicator[l] = true;
                }
                k++;
            }
            
            // check the indicator and shift res_chunk[l] by an appropriate amount
            l = 0;
            while(1){
                if(indicator[l]){   // find the first nonzero entry
                    if(l != 0){
                        evaluator.negacyclic_shift_poly_inplace(res_chunk[l], l * n_test, MemoryPoolHandle().New(false));
                    }
                
                    if(l == (nparallel - 1)){   // there is no more ciphertext to be considered
                        res[j] = res_chunk[l];
                        break;
                    }
                    
                    // check the other remaining terms
                    long l1 = l + 1;
                    while (1) {
                        if(indicator[l1]){
                            evaluator.negacyclic_shift_poly_inplace(res_chunk[l1], l1 * n_test, MemoryPoolHandle().New(false));
                            evaluator.add_inplace(res_chunk[l], res_chunk[l1]);
                        }
                        l1++;
                        if(l1 == nparallel){
                            break;
                        }
                    }
                    res[j] = res_chunk[l];
                    break;
                }
                l++;
            }
            
            // add the (nonzero) intercepts by computing
            // res[j] += sum_{0<=i<ntest} coeff * x^i
            if(scale_and_modt(coeff, model0[j], W0scale, coeff_modulus)){
                evaluator.add_plain_inplace_ckks(res[j], coeff, n_test);
            }
        }
    }
    MT_EXEC_RANGE_END;
}

/*
 Decrypt the input ciphertexts and then decode the resulting vectors.
 
 @param[in] encres[p], the encrypted results
 @param[in] encres_index, the indicator that represents whether the corresponding tag ciphertexts are used or not for performing the imputation at each target variant.
 @param[in] coeff_modulus, the coefficient modulus in the ckks
 @param[in] scale, the scale of the messages
 @param[out] res, the decrypted results
 */
void ckksHEmpute::decrypt_impute(dmat& res, vector<Ciphertext> encres, vector<bool> encres_index,
                                 dvec model0, size_t coeff_modulus, double scale)
{
    long nsnp = res.size();
    long n_test = res[0].size();
    
    MT_EXEC_RANGE(nsnp, first, last);
    for(long j = first; j < last; ++j){
        if(encres_index[j]){
            Plaintext plaintmp;
            decryptor.decrypt_ckks(encres[j], plaintmp);
            uint64_t* plainvec = plaintmp.data();       // hex -> decimal
            
            vector<uint64_t> umsg;
            for(long i = 0; i < n_test; ++i){
                umsg.push_back(plainvec[i]);
            }
            
            vector<double> dmsg;
            decode_Vector(dmsg, umsg, coeff_modulus, scale);
            res[j] = dmsg;
        }
        else {
            vector<double> msg (n_test, model0[j]);
            res[j] = msg;
        }
    }
    MT_EXEC_RANGE_END;
}
