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
#include <sys/resource.h>   // memory usage
#include <iomanip>          // std::setprecision

#include "utils.h"
#include "utils_data.h"
#include "thread.h"

using namespace std;

#define DEBUG false


/*
 @param[in] path, the path for an input file
 @param[in] split_char, the character used for delimiting an input file
 @param[in] header, 1 if the data has a header (so skipping the first row when storing); otherwise 0
 @param[out] res, string-type matrix
 @throws std::invalid_argument if the data is not readable
 */
void SimpleDataFromFile(vector<vector<string>>& res, string path, char split_char, bool header)
{
    ifstream openFile(path.data());
    
    if(openFile.is_open()) {
        string line;
        if(header == true){
            getline(openFile, line);
        }
        /*  read each line */
        while(getline(openFile, line)){
            //vector<string> vecsline;
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            res.push_back(tokens);
            
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
}

/*
 @param[in] path, the path for an input file
 @param[in] split_char, the character used for delimiting an input file
 @param[in] header, 1 if the data has a header (so skipping the first row when storing); otherwise 0
 @param[out] res, double-type matrix
 @throws std::invalid_argument if the data is not readable
 */
void SimpleDataFromFile(vector<vector<double>>& res, string path, char split_char, bool header)
{
    ifstream openFile(path.data());
    
    vector<vector<string>> res_str;
    if(openFile.is_open()) {
        string line;
        
        if(header == true){
            getline(openFile, line);
            //cout << "read header" << endl;
        }
        /*  read each line */
        while(getline(openFile, line)){
            vector<string> vecsline;
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            
            res_str.push_back(tokens);
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
    
    for(long i = 0; i < (long) res_str.size(); ++i){
        vector<double> dtemp;
        for(long j = 0; j < (long) res_str[i].size(); ++j){
            dtemp.push_back(atof(res_str[i][j].c_str()));
        }
        res.push_back(dtemp);
    }
}

/*
 @param[in] path, the path for an input file
 @param[out] res, string-type vector
 @throws std::invalid_argument if the data is not readable
 */
void SimpleDataFromFile(vector<string>& res, string path)
{
    string line;
    
    ifstream myfile (path);
    if (myfile.is_open())
    {
        while (getline (myfile,line))
        {
            res.push_back(line);
        }
        myfile.close();
    }
    else{
        throw invalid_argument("Error: cannot read file");
    }
    
}

/*
 Each row of the data in the path directory is correspoding to one target SNP
 where the entries are the per SNP parameters.
 
 @param[in] path, the file name for trained weights
 @param[out] model_wt, [n_snptarget][2d+1] which is obtained by taking the first (2d+1) lines
 the first column is corresponding to the intercepts of the models
 @param[out] tag_geno_model_coordinates, [n_snptarget][2d] which is obtained by taking the next (2d) lines in the params file
 @param[out] target_geno_model_coordinates, [n_snptarget] where this value is corresponding to the last line
 @throws std::invalid_argument if the data is not readable
*/
void Read_Params(vector<vector<double>> &model_wt, vector<vector<string>> &tag_model_coordinates,
                 vector<string> &target_model_coordinates, string path)
{
    // 1. Read the data and store it as a string-type matrix
    vector<vector<string> > res;
    char split_char = 0x09;
    long nrows = 0;
    ifstream openFile(path.data());
    
    if(openFile.is_open()) {
        string line;
        
        /*  read each line */
        while(getline(openFile, line)){
            istringstream split(line);
            vector<string> tokens;
            for (string each; getline(split, each, split_char); tokens.push_back(each));
            res.push_back(tokens);
            nrows++;
        }
    } else {
        throw invalid_argument("Error: cannot read file");
    }
    
    // 2. Encode each row to weight, tag_model_coordinates, and target_coordinates
    long num_params = res[0].size()/2 ;     // the number of the predictors including the bias term (2*vicinity + 1)
    long dim = num_params - 1;
    long num_params2 = num_params + dim;    // (2v+1) + (2v)
    
    vector<string> original_target_model_coordinates;
    original_target_model_coordinates.resize(nrows);
    
    for(long i = 0; i < nrows; ++i){
        original_target_model_coordinates[i] = (res[i][num_params2]);
    }
    
    // 3.1. read the selected target_geno_model_coordinates
    SimpleDataFromFile(target_model_coordinates, "data/target_geno_model_coordinates.txt");
    long n_snptarget_model = target_model_coordinates.size();   // the actual number of target SNPs (80,882)

   
    // 2. find the starting/ending coordinates
    long n_snptarget_whole = original_target_model_coordinates.size();      // 82993
    long dim1 = dim + 1;
    model_wt.resize(n_snptarget_model, vector<double> (dim1));
    tag_model_coordinates.resize(n_snptarget_model, vector<string> (dim));
    
    string target_starting_var = target_model_coordinates[0]; // the selected coordinate,
    string target_ending_var = target_model_coordinates[n_snptarget_model - 1]; // the selected coordinate,
    
    long jstart = 0;
    while(1){
        if(original_target_model_coordinates[jstart] == target_starting_var){
            break;
        }
        jstart++;
        if(jstart == n_snptarget_whole){
            throw invalid_argument("Error: cannot find the starting");
        }
    }
    
    long jend = n_snptarget_whole - 1;
    while(1){
        if(original_target_model_coordinates[jend] == target_ending_var){
            break;
        }
        jend--;
        if(jend < 0){
            throw invalid_argument("Error: cannot find the ending");
        }
    }
    
    // store the corresponding information of model weights and coordinates
    // jstart <= j <= jend
    long i = 0;
    for(long j = jstart; j <= jend; ++j){
        for(long d = 0; d < dim1; ++d){
            model_wt[i][d] = atof(res[j][d].c_str());
        }
        for(long d = 0; d < dim; ++d){
            tag_model_coordinates[i][d] = res[j][d + dim1];
        }
        //cout << i << ": " << model_wt[i][0] << "," <<  tag_model_coordinates[i][0] << "," << endl;
        //cout << res[j][dim + dim1] << endl;
        i++;
    }
}


/*
 Read the tagSNPs
 (Each row corresponds to a SNP and first 4 columns describe the SNP and remaining columns are the genotypes:
 [Chromosome] [Start] [End] [Name], [Genotype for 1st sample] [Genotype for 2nd sample])
 
 @param[in] path, the file name for an input data
 @param[in] split_char, the character used for delimiting an input file
 @param[out] res_list, [n_snptag]: the list of Starting or Ending
 @param[out] res_array, [n_snptarget][n]: genotype array where n is the number of samples
 @throws std::invalid_argument if the data is not readable
 */
void TagDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array,
                     string path, char split_char)
{
    long list_index = 1;     // use "Start"
    
    ifstream openFile(path.data());
    
    if(openFile.is_open()) {
        string line;
        while(getline(openFile, line)){
            istringstream split(line);
            vector<string> tokens;
            
            long j = 0;
            string each;
            while(getline(split, each, split_char)){
                if(j > 3){
                    tokens.push_back(each);         // real genotype values
                }
                else if(j == list_index){
                    res_list.push_back(each);
                }
                j++;
            }
            res_array.push_back(tokens);
        }
    }
    else {
        throw invalid_argument("Error: cannot find the ending");
    }
}

/*
 Read the targetSNPs
 (Each row corresponds to a SNP and first 4 columns describe the SNP and remaining columns are the genotypes:
 [Chromosome] [Start] [End] [Name], [Genotype for 1st sample] [Genotype for 2nd sample])
 
 @param[in] path, the file name for an input data
 @param[in] array, 1 if we store both list and array; otherwise 0
 @param[in] split_char, the character used for delimiting an input file
 @param[out] res_list, [n_snptag]: the list of Starting or Ending
 @param[out] res_array, [n_snptarget][n]: genotype array where n is the number of samples
 @throws std::invalid_argument if the data is not readable
 */
void TargetDataFromFile(vector<string>& res_list, vector<vector<string>>& res_array,
                        string path, bool array, char split_char)
{
    long list_index = 1;     // use "Start"
    
    
    ifstream openFile(path.data());
    
    // just store the list 
    if(array == 0){
        if(openFile.is_open()) {
            string line;
            while(getline(openFile, line)){
                istringstream split(line);
                vector<string> tokens;
                
                long j = 0;
                string each;
                while(getline(split, each, split_char)){
                    if(j == list_index){
                        res_list.push_back(each);
                        break;
                    }
                    j++;
                }
            }
        }
        else {
            throw invalid_argument("Error: cannot find the ending");
        }
    }
    else{
        if(openFile.is_open()) {
            string line;
            while(getline(openFile, line)){
                istringstream split(line);
                vector<string> tokens;
                
                long j = 0;
                string each;
                while(getline(split, each, split_char)){
                    if(j > 3){
                        tokens.push_back(each);
                    }
                    else if(j == list_index){
                        res_list.push_back(each);
                    }
                    j++;
                }
                
                res_array.push_back(tokens);
            }
        }
        else {
            throw invalid_argument("Error: cannot find the ending");
        }
    }
}

/*
 @param[in] model_data, [n_snptarget][2d+1] which is obtained by taking the first (2d+1) lines in each xxx.params
 @param[in] tag_model_coordinates, [n_snptarget][2d] which is obtained by taking the next (2d) lines in the params file (vector of starting index)
 @param[in] target_model_coordinates, [n_snptarget] which is corresponding to the last line (without ".params")
 @param[in] tag_filename, the filename for tag genotype dataset
 @param[in] target_filename, the filename for target genotype dataset
 @param[in] selection, 1 if "num of target snps = 20K"; 2 if "num of target snps = 40K"; otherwise 3.
 @param[in] acc_test, 1 if "ytest" is needed; otherwise 0

 @param[out] ytest[n_snptag][n_test]
 @param[out] model0[n_snptag], the intercept
 @param[out] model[n_snptag][dim]
 @param[out] tag_data, [n_snptag][n_test] = [16184][1004]
 @param[out] tag_model_starting_index[n_snptag], starting index of model in tag (0 <= x < 16184)
 */
void Read_Genotype(dmat& ytest, dvec& model0, dmat& model,
                   vector<vector<int>>& tag_data, vector<long>& tag_model_starting_index,
                   DATAParam &DATAparam,
                   vector <vector<double>> model_data, vector <vector<string>> tag_model_coordinates, vector <string> target_model_coordinates,
                   string tag_filename, string target_filename, int selection, bool acc_test)
{
    long dim = model_data[0].size() - 1;
    
    // 1. Select a subset of target SNPs (to show the scalability)
    // the real number of target SNPs to be imputed
    long n_snptarget_imputed;
    switch (selection) {
        case 1:
            n_snptarget_imputed = 20000;
            break;
        case 2:
            n_snptarget_imputed = 40000;
            break;
        case 3:
            n_snptarget_imputed = model_data.size();    // 80,882
            break;
    }
    
    // 1.1 First read the raw tagSNPs
    vector<string> tag_data_coordinates;        // [16K]: ending location of tag
    vector<vector<string>> tag_str_data;        // [16K][1004]: columns = ([name] [start] [end] [id], ...)
    TagDataFromFile(tag_data_coordinates, tag_str_data, tag_filename);
    
    // 1.2. Convert the tag data into int-type
    long n_snptag = tag_str_data.size();
    long n_test = tag_str_data[0].size();
    tag_data.resize(n_snptag, vector<int> (n_test));
    
    MT_EXEC_RANGE(n_snptag, first, last);
    for(long j = first; j < last; ++j){
        for(long i = 0; i < (long) tag_str_data[0].size(); ++i){
            tag_data[j][i] = atoi(tag_str_data[j][i].c_str()); // j-th variant of the i-th sample,
        }
    }
    MT_EXEC_RANGE_END;
    
    // 2. Read the targetSNPs
    // if acc_test = true, we read all the data
    // otherwise, just read the second column
    
    vector<string> target_data_coordinates; // [83K]
    vector<vector<string>> target_rawdata;  // [83K][1004]
    
    TargetDataFromFile(target_data_coordinates, target_rawdata, target_filename, acc_test);
    
    long n_snptarget = target_data_coordinates.size();
    cout << "> tag(nsnptag,ntest) = (" << n_snptag << ","  << n_test << ") " << endl;
    cout << "> target(nsnptarget,ntest) = (" << n_snptarget << ","  << n_test << ")" << endl;
    
    // Update parameters
    DATAParam others(dim, n_test, n_snptag, n_snptarget, n_snptarget_imputed, Thread::availableThreads());
    DATAparam = others;

    // 3. store the index of tag/target genotype data index
    long dim1 = dim + 1;
    vector<double> model_onerow (dim, 0.0);
    model0.resize(n_snptarget_imputed);
    model.resize(n_snptarget_imputed, model_onerow);
    tag_model_starting_index.resize(n_snptarget_imputed);
    
    if(acc_test){
        vector<double> ytest_onerow (n_test, 0.0);
        ytest.resize(n_snptarget_imputed, ytest_onerow);
    }
    
    // 4. Divide the whole dataset into nthread-chunks and
    long chunks_size = floor((double) (n_snptarget_imputed)/ DATAparam.n_thread);
    
    MT_EXEC_RANGE(DATAparam.n_thread, first, last); // DATAparam.n_thread
    for(long i = first; i < last; ++i){
        long jstart = i * chunks_size;
        long jend;
        if(i != (DATAparam.n_thread - 1)){
            jend = (i + 1) * chunks_size;
        } else{
            jend = n_snptarget_imputed;
        }
        
        long j_target = 0;
        long j_tag = 0;
        for (long j = jstart; j < jend; ++j){
            // 4.1. read the target model coordinates at j-th entry
            // and find the real current coordinate in the "target_data_coordinates" (e.g., 17084717~50999681)
            string target_model_var = target_model_coordinates[j]; // the coordinate of the current target_model variant
            
            while(1){
                if(target_data_coordinates[j_target] == target_model_var){
                    break;
                }
                j_target++;
                if(j_target == n_snptarget){
                    j_target = 0;
                }
            }
            
            if(acc_test){
                dvec ytest_tmp;     // [ntest]
                for(long k = 0; k < n_test; ++k){
                    ytest_tmp.push_back(atof(target_rawdata[j_target][k].c_str())/2.0);
                }
                ytest[j] = (ytest_tmp);
            }
        
            // 4.2. read the jth row (focusing the first entry) in the tag_model_coordinates
            string tag_starting_var = tag_model_coordinates[j][0];
            
            while(1){
                if(tag_data_coordinates[j_tag] == tag_starting_var)
                    break;
                j_tag++;
                if(j_tag == n_snptag){
                    j_tag = 0;
                }
            }
            tag_model_starting_index[j] = j_tag;
            
            // 4.3. store the model
            dvec model_tmp;
            for(long k = 1; k < dim1; ++k){
                model_tmp.push_back(model_data[j][k]);
            }
            model[j] = (model_tmp);
            model0[j] = (model_data[j][0]);
        }
    }
    MT_EXEC_RANGE_END;
    
    //cout << "ytest = [" << ytest.size() << "][" << ytest[0].size() << "], ";
    //cout << "model0 = (" << model0.size() << "), ";
    cout << "> model = (" << model.size() << ","  << model[0].size() << ")"<< endl;

}
