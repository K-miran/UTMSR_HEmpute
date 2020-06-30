# HEmpute-Test: Secure Imputation via Homomorphic Encryption

HEmpute-Test provides a C++ implementation of secure genotype imputation protocol using the trained linear models from HEmpute-Train. We present two solutions based on homomorphic encryption frameworks - BFV [(https://eprint.iacr.org/2012/144)](https://eprint.iacr.org/2012/144) and CKKS [(https://eprint.iacr.org/2016/421)](https://eprint.iacr.org/2016/421). Our secure imputation protocols are implemented with Microsoft SEAL version 3.4, which includes implementations of BFV and CKKS.

## Table of Contents
- [Data Description](#data-description)
- [Source Code](#source-code)
    - [Requirements](#requirements)
    - [Installation](#installation)
    - [Example Run](#example-run)

## Data Description
### Testing Data
The "/data" directory contains three-types of (compressed) data:
- Tag variant genotypes: tag_testing.txt, tag_testing_AFR.txt, tag_testing_AMR.txt, tag_testing_EUR.txt.
- Target variant genotypes: target_testing.txt, target_testing_AFR.txt, target_testing_AMR.txt, target_testing_EUR.txt.
- target_geno_model_coordinates.txt 

You can also download these data via this [link](https://github.com/K-miran/secure-imputation/tree/master/data).
We excluded the variants at the very end of the chromosome 22 and at the middle of the chromosome (centromere) in the whole target SNPs because we do not have many tag SNPs around those locations. So, the "target_geno_model_coordinates.txt " contains the start coordinates of the target SNPs that were actually used for imputation in our experiment. 

In our protocol, we will input the genotypes in tag_testing.txt to the models and accuracy will be tested using target genotype data in target_testing.txt. 
The genotype files are tab-delimited and each row corresponds to a SNP. First 4 columns describe the SNP and remaining columns are the genotypes:
[Chromosome] [Start] [End] [Name] [Genotype for 1st sample] [Genotype for 2nd sample] ...
Each genotype is coded as 0/1/2. Also, we perform population stratification, so we divide the training and testing samples into 3 super-populations African (AFR), Americans (AMR), and European (EUR). The same directory contains the training/testing tag/target SNP genotypes for these samples.

### Model parameters
The trained parameters are found in the "/params" directory. The files are tab-delimited text files. Each row corresponds to a SNP.

## Source Code

### Requirements
We recommend installing the following: m4, texinfo, homebrew, and cmake. You can easily install them by running the following commands :
```
sudo apt-get install m4
sudo apt-get install texinfo 
mkdir homebrew && curl -L https://github.com/Homebrew/brew/tarball/master | tar xz --strip 1 -C homebrew
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install cmake  
```

### Installation
First, you can install the modified SEAL library the following commands:
```
bash install.sh
```
Then you can install our secure genotype imputation library by running the following commands :
```
cmake . 
make
```

### Example Run
The program will run with the BFV or CKKS homomorphic encryption scheme.
For example, we first make a new folder "res" and run the test program "hefoo"  by running main_hempute.cpp as follows:

```
./hefoo bfv ALL 16 20000 2 est
./hefoo ckks AFR 16 40000 16 label
./hefoo ckks AMR 16 40000 16 microAUC
./hefoo ckks EUR 16 80000 32 macroacc
```
As in the example, the following list of command-line arguments is given after the name of the test program:
- An HE scheme name (e.g. bfv or ckks).
- Data type (e.g. ALL, AFR, AMR, EUR). 
- Number of threads.
- Number of target SNPs (e.g. 20000, 40000, 80000).
- Vicinity size of the imputation linear model (e.g. 2, 4, 8, 16, 24, 32). That is, each variant genotype is modeled using genotypes of variants within variant vicinity of the variant.
- The output format
    - nulll: declare nothing in this scope
    - est: output the predicted estimations
    - label: output the prediction labels for 0,1,2
    - microAUC: output the actual genotypes and  the estimated genotypes in order to calculate the micro-AUCs
    - macroacc: calculate the macro-aggregated accuracies over all variants and non-reference genotypes. 

Then the results are stored in the "/res" directory. 
We also provide the code to measure the micro-AUC. After running the C++ test program with the "micro-AUC" mode, you can obtain the accuracy results by running microAUC_evaluation.py as follows: 
```
./python microAUC_evaluation.py -i inputfile -o outputfile
```
For example, if running the test program "./hefoo ckks ALL 16 80000 16 microAUC", the input list file will be stored as "res/microAUC_ckks_ALL_80000_32.txt ". So, you can run the code by the following command
```
./python microAUC_evaluation.py -i res/microAUC_ckks_ALL_80000_32.txt -o res/microAUC_ckks_ALL_80000_32.png 
```
