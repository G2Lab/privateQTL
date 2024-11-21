#ifndef INPUT_H
#define INPUT_H
 
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <sstream>
#include "utils.h"
using namespace std;

struct GeneData {
    std::string chr;
    uint64_t start;
    uint64_t end;
};

class prepareInput
{
public:
    uint64_t sampleN;
    uint64_t phenoN;
    uint64_t genoN;
    uint64_t ciswindow;
    unordered_map<string, GeneData> genePos;
    vector<uint64_t> snpPos;
    vector<string> snpChr;
    vector<string> snpIDs;
    vector<vector<double>> geno;
    prepareInput(){};
    prepareInput(string& pheno_pos, string& geno_matrix, string& geno_pos, uint64_t window);
    string getCisRange(string geneID, vector<uint64_t>& positions);
    vector<uint64_t> getSNPrange(uint64_t start, uint64_t end, string chr, vector<string>& cisSnpIds);
    int sliceGeno(vector<uint64_t> positions, string& chr, int64_t missing, vector<string>& cisSNPs, vector<vector<double>>& slicedmatrix);
private:
    // Add any private member functions or variables if needed
};

#endif 
