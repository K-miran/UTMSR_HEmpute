#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/resource.h>
#include <iomanip>          // std::setprecision

#include "utils.h"
#include "utils_data.h"
#include "thread.h"
#include "TestHEmpute.h"

// ./hefoo ckks ALL 16 20000 2 mAUC
// ./hefoo ckks AFR 16 20000 16 0 1

using namespace std;
using namespace seal;

int main(int argc, char** argv) {
    
    string scheme = argv[1];            // bfv, ckks
    string data = argv[2];              // e.g. ALL, AFR, AMR, EUR
    long n_thread = atoi(argv[3]);      // e.g. 1, 2, ..., 24
    long n_target = atoi(argv[4]);      // e.g. 20000, 40000, 80000
    long vicinity = atoi(argv[5]);      // e.g. 2, 4, 8, 16, 24, 32
    string acc_test = argv[6];          // label, est, microAUC, macroacc

    long dim = 2 * vicinity;
    
    cout << "+------------------------------------+" << endl;
    cout << "|              PARAMETERS            |" << endl;
    cout << "+------------------------------------+" << endl;
    
    cout << "> Scheme: " << scheme << ", nthread: " << n_thread << ", ntargets: " << n_target << ", npredictors: " << dim << endl;
    
    int selection = 0;
    switch (n_target) {
        case 20000:
            selection = 1;
            break;
        case 40000:
            selection = 2;
            break;
        case 80000:
            selection = 3;
            break;
    }
    
    Thread::initThreadPool(n_thread);

    // model parameters
    vector<vector<double>> model_data;                // [n_snptarget][2d+1] which is obtained by taking the first (2d+1) lines in each xxx.params
    vector<vector<string>> tag_model_coordinates;     // [n_snptarget][2d] which is obtained by taking the next (2d) lines in the params file
    vector<string> target_model_coordinates;          // [n_snptarget] where this value is corresponding to the last line (without ".params")
    string model_path;
    
    dmat ytest;     // the actual genotypes
    dmat HE_ypred;  // estimated genotypes
    
    struct rusage usage;
    chrono::high_resolution_clock::time_point time_start, time_end;
    
    if(data == "ALL"){
        cout << "+------------------------------------+" << endl;
        cout << "|         0.1. Read the model        |" << endl;
        cout << "+------------------------------------+" << endl;

        time_start = chrono::high_resolution_clock::now();
        
        // 0.1 Read the trained model parameters from the docker
        if(vicinity == 32){
            model_path = "params/SAVED_PARAMS_1_51000000_32.params";
        } else{
            model_path = "params/LMSE_models_" + to_string(vicinity) + ".params";
        }
            
        Read_Params(model_data, tag_model_coordinates, target_model_coordinates, model_path);
            
        time_end = chrono::high_resolution_clock::now();
        auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
        cout << "Select " << model_data.size() << " models with " << model_data[0].size() - 1 << "-predictors (seconds) : " << setprecision(4) << time_diff.count()/(1000.0)  << endl;
        getrusage(RUSAGE_SELF, &usage);
        cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

        cout << "+------------------------------------+" << endl;
        cout << "|     0.2. Read the tag & target     |" << endl;
        cout << "+------------------------------------+" << endl;

        time_start = chrono::high_resolution_clock::now();
        string tag_filename = "data/tag_testing.txt";           // [ptag][n+4] = [16184][1004 + 4]
        string target_filename = "data/target_testing.txt";     // [ptarget][n+4] = [83072][1004 + 4]

        dvec model0;    // the intercepts of each model
        dmat model;     // the trained model parameters

        DATAParam DATAparam;
        vector<vector<int>> tag_geno_data;      // [p][n], the tag genotype data
        vector<long> tag_model_starting_index;   // starting index of model in tag
        bool is_ytest_needed;
        if((acc_test == "microAUC")||(acc_test == "macroacc")){
            is_ytest_needed = true;
        } else{
            is_ytest_needed = false;
        }
        
        time_start = chrono::high_resolution_clock::now();

        Read_Genotype(ytest, model0, model,
                      tag_geno_data, tag_model_starting_index, DATAparam,
                      model_data, tag_model_coordinates, target_model_coordinates,
                      tag_filename, target_filename, selection, is_ytest_needed);

        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
        cout << "Curate the model & dataset (seconds) : " << time_diff.count()/(1000.0) << endl;
        getrusage(RUSAGE_SELF, &usage);
        cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;

        cout << "+------------------------------------+" << endl;
        cout << "|          1. HE-Prediction          |" << endl;
        cout << "+------------------------------------+" << endl;
        TestHEmpute testHEmpute(DATAparam);
        HE_ypred.resize(DATAparam.n_snptarget_model, vector<double> (DATAparam.n_test)); // [p][n]
        
        if(scheme == "bfv"){
            testHEmpute.bfv_HEmpute(HE_ypred, model0, model, tag_geno_data, tag_model_starting_index);
        }
        else if(scheme == "ckks"){
            testHEmpute.ckks_HEmpute(HE_ypred, model0, model, tag_geno_data, tag_model_starting_index);
        }

        //print_data(ytest, HE_ypred, predmAUC_filename);
    }
    else {      // This is for the population stratified imputation
        cout << "+------------------------------------+" << endl;
        cout << "|         0.1. Read the model        |" << endl;
        cout << "+------------------------------------+" << endl;
        
        time_start = chrono::high_resolution_clock::now();
        
        model_path = "params/SAVED_PARAMS_1_51000000_32_" + data + ".params";
        
        Read_Params(model_data, tag_model_coordinates, target_model_coordinates, model_path);
        
        time_end = chrono::high_resolution_clock::now();
        auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
        cout << "Select " << model_data.size() << " many models with " << model_data[0].size() - 1 << "-predictors (seconds) : " << time_diff.count()/(1000.0)  << endl;
        getrusage(RUSAGE_SELF, &usage);
        cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
        
        cout << "+------------------------------------+" << endl;
        cout << "|     0.2. Read the tag & target     |" << endl;
        cout << "+------------------------------------+" << endl;
        
        string tag_filename = "data/tag_testing_" + data + ".txt";          // [ptag][n+4] = [16184][1004 + 4]
        string target_filename = "data/target_testing_" + data + ".txt";    // [ptarget][n+4] = [83072][1004 + 4]
        
        bool is_ytest_needed;
        if((acc_test == "microAUC")||(acc_test == "macroacc")){
            is_ytest_needed = true;
        } else{
            is_ytest_needed = false;
        }
        
        time_start = chrono::high_resolution_clock::now();
        
        dvec model0;    // the intercepts of each model
        dmat model;     // the trained model parameters
        
        DATAParam DATAparam;
        vector<vector<int>> tag_geno_data;      // [p][n], the tag genotype data
        vector<long> tag_model_starting_index;   // starting index of model in tag
        
        time_start = chrono::high_resolution_clock::now();
        
        Read_Genotype(ytest, model0, model,
                      tag_geno_data, tag_model_starting_index, DATAparam,
                      model_data, tag_model_coordinates, target_model_coordinates,
                      tag_filename, target_filename, selection, is_ytest_needed);
        
        time_end = chrono::high_resolution_clock::now();
        time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
        cout << "Curate the model & dataset (seconds) : " << time_diff.count()/(1000.0) << endl;
        getrusage(RUSAGE_SELF, &usage);
        cout<< "RAM Usage (GB): "  << (double) usage.ru_maxrss/(DATAParam::memoryscale) << endl;
        
        cout << "+------------------------------------+" << endl;
        cout << "|          1. HE-Prediction          |" << endl;
        cout << "+------------------------------------+" << endl;
        
        TestPOPHEmpute testHEmpute(DATAparam);
        HE_ypred.resize(DATAparam.n_snptarget_model, vector<double> (DATAparam.n_test)); // [p][n]
        bool parallel = true;
        
        if(scheme == "bfv"){
            testHEmpute.bfv_HEmpute(HE_ypred, model0, model, tag_geno_data, tag_model_starting_index, parallel);
        }
        else if(scheme == "ckks"){
            testHEmpute.ckks_HEmpute(HE_ypred, model0, model, tag_geno_data, tag_model_starting_index, parallel);
        }
    }
    
    cout << "+------------------------------------+" << endl;
    cout << "|             2. Accuracy            |" << endl;
    cout << "+------------------------------------+" << endl;
    
    time_start = chrono::high_resolution_clock::now();
    
    if(acc_test == "label"){
        string filename = "res/" + scheme  + "_" + data + "_" + to_string(n_target) + "_" + to_string(dim) + ".txt";
        ofstream outf(filename);
        outf.close();
        print_labels(filename, HE_ypred, target_model_coordinates);
    }
    else if(acc_test == "est"){
        string filename = "res/estimates_" + scheme  + "_" + data + "_" + to_string(n_target) + "_" + to_string(dim) + ".txt";
        ofstream outf(filename);
        outf.close();
        print_estimates(filename, HE_ypred, target_model_coordinates);
    }
    else if(acc_test == "microAUC"){
        string filename = "res/microAUC_" + scheme  + "_" + data + "_" + to_string(n_target) + "_" + to_string(dim) + ".txt";
        print_data(ytest, HE_ypred, filename);
    }
    else if(acc_test == "macroacc"){
        double allmacro = test_allvariant_macro_acc(ytest, HE_ypred);
        cout << "> All Variant macro accuracy = " << allmacro << endl;
        double nonrefmacro = test_nonrefgenotype_macro_acc(ytest, HE_ypred);
        cout << "> Non-ref Genotype macro accuracy = " << nonrefmacro << endl;
    }
    
    time_end = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::milliseconds>(time_end - time_start);
    cout << "Print the predicted results (seconds) : " << time_diff.count()/(1000.0)  << endl;
    
    return 0;
}
