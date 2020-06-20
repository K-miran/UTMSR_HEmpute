/*
 * @file       TestHEmpute.h, header file
 * @brief      defining functions for genotype imputation
 *
 * @author     Miran Kim
 * @date       June. 4, 2019
 * @copyright  GNU Pub License
 */


#include <vector>
#include "utils.h"
#include "seal/seal.h"

using namespace std;
using namespace seal;

class TestHEmpute
{
public:
    
    DATAParam& DATAparam;
    
    TestHEmpute(DATAParam& DATAparam): DATAparam(DATAparam) {}

    void bfv_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data, vector<long> tag_index);
    
    void ckks_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data, vector<long> tag_index);
    
};

class TestPOPHEmpute
{
public:
    
    DATAParam& DATAparam;
    
    TestPOPHEmpute(DATAParam& DATAparam): DATAparam(DATAparam) {}
    
    void bfv_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data, vector<long> tag_index, bool parallel);
    
    void ckks_HEmpute(dmat& ypred, dvec model0, dmat model, vector<vector<int>> tag_geno_data, vector<long> tag_index, bool parallel);
};
