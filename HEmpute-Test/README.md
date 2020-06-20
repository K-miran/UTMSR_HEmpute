# HEmpute-Evaluation

### Table of Contents
- [Data Description](#data-description)
- [Source Code](#source-code)
    - [Requirements](#requirements)
    - [Building the code](#building-the-code)
    - [Running the code](#running-the-code)

## Data Description
### Testing Data
The "/data" directory contains three-types of (compressed) data:
- Tag SNP genotypes: tag_testing.txt, tag_testing_AFR.txt, tag_testing_AMR.txt, tag_testing_EUR.txt.
- Target SNP genotypes: target_testing.txt, target_testing_AFR.txt, target_testing_AMR.txt, target_testing_EUR.txt.
- target_geno_model_coordinates.txt 

You can also download these data via this [link](https://github.com/K-miran/secure-imputation/tree/master/data).
We excluded the variants at the very end of the chromosome 22 and at the middle of the chromosome (centromere) in the whole target SNPs because we do not have many tag SNPs around those locations. So, the "target_geno_model_coordinates.txt " contains the start coordinates of the target SNPs that were actually used for imputation in our experiment. 
In addition, we provide "target_geno_model_coordinates_ending.txt " which contains the end coordinates of the actual target SNPs. 

In our protocol, we will input the genotypes in tag_testing.txt to the models and accuracy will be tested using target_testing.txt genotype data. 
The genotype files are tab-delimited and each row corresponds to a SNP. First 4 columns describe the SNP and remaining columns are the genotypes:
- [Chromosome] [Start] [End] [Name] [Genotype for 1st sample] [Genotype for 2nd sample] ...

Each genotype is coded as 0/1/2. Also, we performed population stratification, so we divided the training and testing samples into 3 super-populations African (AFR), Americans (AMR), and European (EUR). The same directory contains the training/testing tag/target SNP genotypes for these samples.

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

### Building the code
First, you can install SEAL library the following commands:
```
bash install.sh
```
Then you can install our secure genotype imputation library by running the following commands :
```
cmake . 
make
```

### Running the code
The program will run with the BFV or CKKS homomorphic encryption scheme.
For example, we run the test program "hefoo"  by running main_heimpute.cpp as follows:

```
./hefoo bfv ALL 16 20000 2 est
./hefoo ckks ALL 16 40000 16 label
./hefoo ckks EUR 16 80000 32 mAUC
```
As in the example, the following list of command-line arguments is given after the name of the test program:
- An HE scheme name (e.g. bfv or ckks).
- Data type (e.g. ALL, AFR, AMR, EUR). 
- Number of threads.
- Number of target SNPs (e.g. 20000,40000,80000).
- Number of vicinities of the imputation linear model (e.g. 2,4,8,16,24,32).
- The output format: est if you want to output the predicted estimations; label if you want to output the prediction labels for 0,1,2; mAUC if you want to calculate the micro-AUCs of the actual genotypes and the estimated genotypes.
