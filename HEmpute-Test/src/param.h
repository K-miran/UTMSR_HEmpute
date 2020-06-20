/*
 * @file       param.h, header file
 * @brief      defining functions for parameter
 *
 * @author     Miran Kim
 * @date       June. 3, 2019
 * @copyright  GNU Pub License
 */

#include <iostream>
#include <vector>
#include <string>
#include "seal/seal.h"

using namespace std;
using namespace seal;

class DATAParam {
public:
    long dim;
    long n_test;
    long n_snptag;
    long n_snptarget;
    long n_snptarget_model; // in the whole testing, n_snptarget = n_snptarget_model
    long n_thread;
    
    static long memoryscale;
    
    // default constructor
    DATAParam() {}
    
    // parameterized constructor
    DATAParam(long dim, long n_test, long n_snptag, long n_snptarget, long n_snptarget_model, long n_thread): dim(dim), n_test(n_test), n_snptag(n_snptag), n_snptarget(n_snptarget), n_snptarget_model(n_snptarget_model), n_thread(n_thread) {}
    
    void operator=(const DATAParam &other){
        dim = other.dim;
        n_test = other.n_test;
        n_snptag = other.n_snptag;
        n_snptarget = other.n_snptarget;
        n_snptarget_model = other.n_snptarget_model;
        n_thread = other.n_thread;
    }
};



class HEParam {
public:
    long logq;
    long nslots;
    long nparallel;          // maximal number of testing rows, which are performed in parallel
    long nctxts;             // total number of generated ciphertexts for each variant
    long nparallel_last;     // the number of the last rows
    
    // default constructor
    HEParam() {}
    
    HEParam(long logq, long nslots, long nparallel, long nctxts, long nparallel_last): logq(logq), nslots(nslots), nparallel(nparallel), nctxts(nctxts), nparallel_last(nparallel_last) {}
    
    // parameterized constructor
    void operator=(const HEParam &other){
        logq = other.logq;
        nslots = other.nslots;
        nparallel = other.nparallel;
        nctxts = other.nctxts;
        nparallel_last = other.nparallel_last;
    }
};

void getHEParam(HEParam &HEparam, std::shared_ptr<seal::SEALContext> context, long dim, long n_test, long n_snp);


class POPDATAParam {
public:
    long dim;
    vector<long> n_test;
    long n_snptag;
    long n_snptarget;
    long n_snptarget_model; // in the whole testing, n_snptarget = n_snptarget_model
    long n_thread;
    
    static long memoryscale;
    // default constructor
    POPDATAParam() {}
    
    // parameterized constructor
    POPDATAParam(long dim, vector<long> n_test, long n_snptag, long n_snptarget, long n_snptarget_model, long n_thread): dim(dim), n_test(n_test), n_snptag(n_snptag), n_snptarget(n_snptarget), n_snptarget_model(n_snptarget_model), n_thread(n_thread) {}
    
    void operator=(const POPDATAParam &other){
        dim = other.dim;
        for(long i = 0; i < (long) other.n_test.size(); ++i){
            n_test.push_back(other.n_test[i]);
        }
        n_snptag = other.n_snptag;
        n_snptarget = other.n_snptarget;
        n_snptarget_model = other.n_snptarget_model;
        n_thread = other.n_thread;
    }
};

class POPHEParam {
public:
    long logq;
    long nslots;
    vector<long> nparallel;          // maximal number of testing rows, which are performed in parallel
    vector<long> nctxts;             // total number of generated ciphertexts for each variant
    vector<long> nparallel_last;     // the number of the last rows
    
    // default constructor
    POPHEParam() {}
    
    POPHEParam(long logq, long nslots, vector<long> nparallel, vector<long> nctxts, vector<long> nparallel_last): logq(logq), nslots(nslots), nparallel(nparallel), nctxts(nctxts), nparallel_last(nparallel_last) {}
    
    // parameterized constructor
    void operator=(const POPHEParam &other){
        logq = other.logq;
        nslots = other.nslots;
        nparallel = other.nparallel;
        nctxts = other.nctxts;
        nparallel_last = other.nparallel_last;
    }
};

void getHEParam(POPHEParam &HEparam, std::shared_ptr<seal::SEALContext> context, vector<long> n_test, long n_snp);
