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
privateQTL-I assumes phenotype is publicly available and therefore preprocessing is completed in plaintext. It takes in genotype that has been locally projected onto reference PCs and residualized, and fully preprocessed phenotype from each data owner. It is run on a per-gene basis, and can run two gene ranges in parallel (start-middle, middle-end). It takes in the file path for genotype and phenotype, as well as position matrices for indexing. Please run as the following.
```sh
./eqtl_mapping [start_gene_index] [middle_gene_index] [end_gene_index] [num_permutations] [pheno_file_path] [geno_file_path] [pheno_pos] [geno_pos][cis_output_prefix] [nominal_output_prefix]
```

### privateQTL-II: private genotype, private phenotype
privateQTL-II assumes phenotypes are also private, and requires additional MPC phenotype preprocessing. It takes in two ranges of gene indices to run in parallel, pre-computed zscore file with total number of samples across sites, normalization method, output file path, and number of samples in each site. Please run as the following.
```sh
./preprocessing [start_gene_index] [middle_gene_index] [end_gene_index] [pheno_input][zscores_file] [normalization] [output_path] [siteA_n] [siteB_n] [siteC_n]
```
Once the preprocessing has finished, eQTL mapping can be run in the same way as privateQTL-I. 

## Test running from example dataset from GEUVADIS
For convenience, we have provided a toy dataset from the GEUVADIS dataset consisting of 30 samples. Please feel free to test both privateQTL versions. Please note that smaller samples sizes means lower accuracy for privateQTL-II MPC preprocessing. 
Please unzip ```toydata.tar.gz``` file. Inside, you will find:
- toy_GEUVADIS_preprocessed_geno.tsv : genotype data, fully processed via projection onto reference panel.
- toy_GEUVADIS_deseq2_pheno.tsv : phenotype data, fully processed with deseq2 normalization and aggregated PCA.
- toy_pheno_position.tsv : phenotype position indexing.
- toy_geno_position.tsv : genotype position indexing.
- toy_GEUVADIS_raw_reads_pheno.tsv : raw reads, used as input for privateQTL-II preprocessing.
- zscores.txt : pre-computed zscores for 30 samples. Used as input for privateQTL-II preprocessing. 

**For privateQTL-II preprocessing**
The following command will produce ```mpc_preprocessed.tsv``` file, which has been deseq2 normalized, locally corrected with PCA, and inverse normal transformed. 
```sh
./preprocessing 0 8120 16241 \
./toy_GEUVADIS_raw_reads_pheno.tsv \
./zscores.txt \
deseq2 \
./mpc_preprocessed \
10 10 10
```

**For privateQTL-I and II eQTL mapping**
The following command will do the mapping. For privateQTL-I, we provide a phenotype that has been preprocessed in plaintext. For privateQTL-II, please use the output of the MPC preprocessing. The example command is for privateQTL-I, for gene 0 to 10.
```sh
./eqtl_mapping 0 5 10 1000 \
./toy_GEUVADIS_deseq2_pheno.tsv \
./toy_GEUVADIS_preprocessed_geno.tsv \
./toy_pheno_position.tsv \
./toy_geno_position.tsv \
./test \
./testn
```
