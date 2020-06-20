/*
 * @file       param.cpp, cpp file
 * @brief      defining functions for paramete
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
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <memory>
#include <limits>

#include "param.h"

using namespace std;

long DATAParam::memoryscale = (1 << 20);            // (GB) 2^20: linux, 2^30: mac
long POPDATAParam::memoryscale = (1 << 20);
//size_t POPDATAParam::n_thread = 4;

void getHEParam(HEParam &HEparam, std::shared_ptr<seal::SEALContext> context, long dim, long n_test, long nsnp)
{
    auto &context_data = *context->key_context_data();
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
        case seal::scheme_type::BFV:
            scheme_name = "BFV";
            break;
        
        case seal::scheme_type::CKKS:
            scheme_name = "CKKS";
            break;
        
        default:
            throw std::invalid_argument("unsupported scheme");
    }
    
    long logq = 0;
    long nslots = 0;
    long nparallel = 0;
    long nctxts = 0 ;
    long nparallel_last = 0;
    
    if(scheme_name == "BFV"){
        logq = context_data.total_coeff_modulus_bit_count();
        nslots = context_data.parms().poly_modulus_degree();
    }
    else if(scheme_name == "CKKS"){
        logq = context_data.parms().coeff_modulus()[0].value();
        nslots = context_data.parms().poly_modulus_degree();
    }
    
    nparallel = floor((double)(nslots)/(double)(n_test));   // maximal number of testing rows, which are performed in parallel
    nctxts = ceil((double)(nsnp)/(double)(nparallel));      // total number of generated ciphertexts for each variant
    nparallel_last = nsnp - nparallel * (nctxts - 1);     // the number of the last rows
    
    HEParam other(logq, nslots, nparallel, nctxts, nparallel_last);
    HEparam = other;
}

void getHEParam(POPHEParam &HEparam, std::shared_ptr<seal::SEALContext> context, vector<long> n_test, long nsnp)
{
    auto &context_data = *context->key_context_data();
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
        case seal::scheme_type::BFV:
            scheme_name = "BFV";
            break;
            
        case seal::scheme_type::CKKS:
            scheme_name = "CKKS";
            break;
            
        default:
            throw std::invalid_argument("unsupported scheme");
    }
    
    long logq = 0;
    long nslots = 0;
    vector<long> nparallel;
    vector<long> nctxts;
    vector<long> nparallel_last;
    
    if(scheme_name == "BFV"){
        logq = context_data.total_coeff_modulus_bit_count();
        nslots = context_data.parms().poly_modulus_degree();
    }
    else if(scheme_name == "CKKS"){
        logq = context_data.parms().coeff_modulus()[0].value();
        nslots = context_data.parms().poly_modulus_degree();
    }
    
    // we will encrypt "nparallel" columns into a single ciphertext 
    for(long r = 0; r < (long) n_test.size(); ++r){
        nparallel.push_back(floor((double)(nslots)/(double)(n_test[r])));   // maximal number of testing rows, which are performed in parallel
        nctxts.push_back(ceil((double)(nsnp)/(double)(nparallel[r])));      // total number of generated ciphertexts for each variant
        nparallel_last.push_back(nsnp - nparallel[r] * (nctxts[r] - 1));  // the number of the last rows
    }
    
    POPHEParam other(logq, nslots, nparallel, nctxts, nparallel_last);
    HEparam = other;
}
