#ifndef PRIVATEQTL_DATACLIENT_H
#define PRIVATEQTL_DATACLIENT_H
#include "input.h"

void dataclient_mapping(string pheno_path, string geno_path, string pheno_pos, string geno_pos, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend, int permut);
// void dataclient1(string pheno_path, string geno_path, string pheno_pos, string geno_pos, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend, int permut);
void dataclient_preprocessing(string norm_method, string pheno_input, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,string zscorefile, vector<vector<double>>& resultVec, vector<string>& gene_string);
// void client_initialize(int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3);
class dataclient
{
public:
    IOService ios;
    uint64_t p;
    dataclient(){};
    void initialize(int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3);
    void do_recv_string(vector<string>& stringvec);
    void do_send_string(vector<string>& tosend);
    void find_common(string phenofile);
    prepareInput load_genotype(string pheno_pos, string geno_path, string geno_pos, int cis_window);
    void mapping(string pheno_path, string geno_path, string pheno_pos, string geno_pos, int rowstart, int rowend, int permut);
    void phenopreprocess(string norm_method, string pheno_input, int rowstart, int rowend, string zscorefile, vector<vector<double>>& resultVec, vector<string>& gene_string);
    void close();
private:
    Channel owner_p1;
    Channel owner_p2;
    Channel owner_p3;
    Channel p1_owner;
    Channel p2_owner;
    Channel p3_owner;
    prepareInput myinput;
    vector<vector<double>> pheno;
    vector<string> genes;
};
#endif