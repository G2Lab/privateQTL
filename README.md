# privateQTL: Secure and Federated Quantitative Trait Loci Mapping with privateQTL
## Overview
privateQTL is a novel framework for secure and federated cis-eQTL mapping across multiple institutions by leveraging Multiparty Computation (MPC). 

## Installation 
privateQTL uses [cryptoTools](https://github.com/ladnir/cryptoTools) v1.10.1, [NTL](https://libntl.org/doc/tour-unix.html), and [gmp](https://gmplib.org/manual/Installing-GMP) . Install them and make sure to add their directory to the CMakeLists.txt file as the following.
```sh
set(cryptoTools_DIR "/path/to/cryptoTools/cmake")
```
The following dependencies should also be accessible:
- cmake minimum v3.18
- gcc v11.2.0
- Eigen v3.4.0
- gperftools v2.5.91
- gsl v2.6
- openblas v0.2.19
- openmp v9.0.1

To build privateQTL, please run the following.
```sh
mkdir build 
cd build
cmake ..
make
```

## Running privateQTL
privateQTL-I and II shares eQTL mapping code, that takes in secretly shared genotype and phenotype and performs matrix multiplication. privateQTL-II requires additional phenotype preprocessing code in MPC. 
### privateQTL-I: private genotype, public phenotype
privateQTL-I assumes phenotype is publicly available and therefore preprocessing is completed in plaintext. It takes in genotype that has been locally projected onto reference PCs and residualized, and fully preprocessed phenotype from each data owner. It is run on a per-gene basis, and can run two gene ranges in parallel (start-middle, middle-end). Please run as the following.
```sh
./eqtl_mapping [start_gene_index] [middle_gene_index] [end_gene_index] [num_permutations] [pheno_file_path] [data_split_set] [cis_output_prefix] [nominal_output_prefix]
```

### privateQTL-II: private genotype, private phenotype
privateQTL-II assumes phenotypes are also private, and requires additional MPC phenotype preprocessing. It takes in two ranges of gene indices to run in parallel, pre-computed zscore file with total number of samples across sites, normalization method, and split set. Please run as the following.
```sh
./preprocessing [start_gene_index] [middle_gene_index] [end_gene_index] [zscores_file] [normalization] [split_set]
```
Once the preprocessing has finished, eQTL mapping can be run in the same way as privateQTL-I.