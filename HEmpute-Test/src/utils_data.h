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
                 vector <string> &target_model_coordinates, string path);

/*------------------------
 Genotype data
 ------------------------*/

void TagDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array, string path,
                     const char split_char = 0x09);

void TargetDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array, string path,
                        const bool array = false, const char split_char = 0x09);

void Read_Genotype(dmat& ytest, dvec& model0, dmat& model,
                   vector<vector<int>>& tag_geno_data, vector<long>& tag_model_starting_index, DATAParam &DATAparam,
                   vector <vector<double>> model_data, vector <vector<string>> tag_model_coordinates, vector <string> target_model_coordinates,
                   string tag_filename, string target_filename, int selection, const bool acc_test = false);













#endif /* DATABASE_H_ */
