/*
 * @file       utils.h, header file
 * @brief      defining functions for reading and storing data
 *
 * @author     Miran Kim
 * @date       June. 3, 2019
 * @copyright  GNU Pub License
 */

#ifndef UTILSDATA_H_
#define UTILSDATA_H_
 
#include <iostream>
#include <vector>
#include <string>
#include "param.h"
#include "utils.h"

using namespace std;

/*------------------------
 Simple IO
 ------------------------*/

void SimpleDataFromFile(vector<string>& res, string filename);
void SimpleDataFromFile(vector<vector<string> >& sline, string path, const char split_char = 0x09, const bool header = false);
void SimpleDataFromFile(vector<vector<double>>& res, string path, const char split_char = 0x09, const bool header = false);

/*------------------------
 Params data
 ------------------------*/

void Read_Params(vector <vector<double>> &model_wt, vector <vector<string>> &tag_model_coordinates,
                 vector <string> &target_model_coordinates, string path, const bool ending = true);

/*------------------------
 Genotype data
 ------------------------*/

void TagDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array, string path,
                     const bool ending = true, const char split_char = 0x09);

void TargetDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array, string path,
                        const bool ending = true, const bool array = false, const char split_char = 0x09);

void Read_Genotype(dmat& ytest, dvec& model0, dmat& model,
                   vector<vector<int>>& tag_geno_data, vector<long>& tag_model_starting_index, DATAParam &DATAparam,
                   vector <vector<double>> model_data, vector <vector<string>> tag_model_coordinates, vector <string> target_model_coordinates,
                   string tag_filename, string target_filename, int selection, const bool ending = true, const bool acc_test = false);













void SelectWeight(vector<vector<double>> &new_model_wt, vector<vector<string>> &new_tag_model_coordinates,
                  vector<string> &new_target_model_coordinates,
                  vector<vector<double>> model_wt, vector<vector<string>> tag_model_coordinates,
                  vector<string> target_model_coordinates, const bool ending = true);



void Read_Genotype(vector<vector<vector<double>>>& ytest, vector<vector<double>>& model0, vector<vector<vector<double>>>& model,
                   vector<vector<vector<int>>> &tag_geno_data, vector<vector<long>>& tag_model_starting_index,
                   vector<vector<string>> &target_geno_id, POPDATAParam &POPDATAparam,
                   vector<vector<vector<double>>> model_data,
                   vector<vector<vector<string>>> tag_model_coordinates, vector<vector<string>> target_model_coordinates,
                   int selection, vector<string> tag_filename, vector<string> target_filename);




void CurateData(dmat& ytest, dvec& model0, dmat& model,
                vector<vector<int>> &tag_geno_data, vector<long>& tag_index, vector<string> &target_geno_id,
                vector<vector<string>> tag_geno_rawdata, vector<vector<string>> target_geno_data, vector <vector<double>> model_data,
                vector <vector<string>> tag_geno_model_coordinates, vector <string> target_geno_model_coordinates,
                long dim, long n_test, long n_snptag, long n_snptarget, long n_snptarget_model, size_t n_thread);


class PLAINIMPUTE{
public:
    DATAParam &DATAparam;
    
    PLAINIMPUTE(DATAParam &DATAparam): DATAparam(DATAparam) {}
    
    void CurateData(dten& Xtest, dmat& ytest, dmat& model0, dmat& model,
                    vector<vector<string>> tag_geno_data, vector<vector<string>> target_geno_data, vector <vector<double>> model_data,
                    vector <vector<string>> tag_geno_model_coordinates, vector <string> target_geno_model_coordinates);

    void CurateData(vector<long>& tag_index, dmat& ytest, dvec& model0, dmat& model,
                    vector<vector<string>> tag_geno_data, vector<vector<string>> target_geno_data, vector <vector<double>> model_data,
                    vector <vector<string>> tag_geno_model_coordinates, vector <string> target_geno_model_coordinates, long nsnp_trial);
    
    void test_impute_linreg(dmat& ypred, dten Xtest, dmat model0, dmat model);
};


#endif /* DATABASE_H_ */
