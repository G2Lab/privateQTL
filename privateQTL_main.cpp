
#include "privateQTL_cp.h"
#include "privateQTL_dataclient.h"
#include <cmath>
#include "utils.h"
#include <chrono>
#include "input.h"
#include <sys/resource.h>
struct rusage r_usage;

void setThreadAffinity(std::thread& thread, int coreId)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(coreId, &cpuset);
    pthread_setaffinity_np(thread.native_handle(), sizeof(cpu_set_t), &cpuset);
}

bool isPortOpen(int port) {
    int sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) {
        perror("socket");
        return false;
    }

    sockaddr_in addr{};
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    addr.sin_addr.s_addr = htonl(INADDR_ANY);

    int result = bind(sockfd, (struct sockaddr*)&addr, sizeof(addr));
    if (result < 0) {
        close(sockfd);
        return false;
    }

    close(sockfd);
    return true;
}
void startMPCset1(string pheno_path, string geno_path, string pheno_pos, string geno_pos,vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, int CPU_core, int permut, Logger& cislogger, Logger& nominalLogger)
{
    vector<thread> threads;
    thread dataClientThread(dataclient_mapping, pheno_path, geno_path, pheno_pos, geno_pos, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
    availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, permut);
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    // int startidx1 = startidx;
    thread runMPC1([=, &cislogger, &nominalLogger]() {
        cp_mapping(0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));

    thread runMPC2([=, &cislogger,&nominalLogger]() {
        cp_mapping(1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow,permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=, &cislogger,&nominalLogger]() {
        cp_mapping(2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow, permut,cislogger,nominalLogger);
    });
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
}
void startMPCset2(string norm_method, string pheno_input, vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, string zscorefile, int CPU_core, vector<vector<double>>& resultVec, vector<string>& gene_string)
{
    vector<thread> threads;
    thread dataClientThread([=,  &resultVec, &gene_string]() {
        dataclient_preprocessing(norm_method, pheno_input, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
        availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, zscorefile, resultVec, gene_string);
    });
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    thread runMPC1([=]() {
        cp_preprocess(norm_method, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));
    thread runMPC2([=]() {
        cp_preprocess(norm_method, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=]() {
        cp_preprocess(norm_method, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
}
// privateQTL : eQTL mapping (shared by privateQTL-I and privateQTL-II)
int main1(int argc, char* argv[])
{
    if (argc < 8)
    {
        cout << "Please provide at least four arg: low, mid, high, permutation, norm method, cislog, nominallog.\n";
        return 1;
    }
    int lo_row = stoi(argv[1]);
    int mid_row = stoi(argv[2]);
    int hi_row = stoi(argv[3]);
    int permut = stoi(argv[4]);
    // string zscorefile = argv[5]; // don't need zscores for scenario1
    string pheno_path = argv[5];
    string geno_path = argv[6];
    string pheno_pos= argv[7];
    string geno_pos = argv[8];
    string cis_log = argv[9];
    string nominal_log = argv[10];
    int startingPort = 12300;
    int assignedPorts=0;
    vector<int> openPorts;
    while (assignedPorts < 24)
    {
        int port = startingPort;
        if(isPortOpen(port))
        {
            openPorts.push_back(port);
            assignedPorts++;
        }
        startingPort++;
    }
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    int numCores = thread::hardware_concurrency();
    vector<thread> threads;
    // string nominal = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/" + split_set + "/privateQTL_scenario1_"+split_set+"_"+nominal_log);
    // string cis = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/output/" + split_set + "/privateQTL_scenario1_"+split_set+"_"+cis_log);
    Logger nominallogger(nominal_log+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv"),  cislogger(cis_log+"_"+to_string(lo_row)+"_"+to_string(mid_row)+".tsv");
    Logger nominallogger2(nominal_log+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv"), cislogger2(cis_log+"_"+to_string(mid_row)+"_"+to_string(hi_row)+".tsv");
    nominallogger.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    nominallogger2.log(string("phenID\tvarID\tvarIdx\tdof\tr_nom\tr2_nom\ttstat\tpval\tslope\tslope_se"));
    cislogger2.log(string("phenID\tvarID\tvarIdx\tbeta_shape1\tbeta_shape2\ttrue_dof\tpval_true_df\tr_nom\tr2_nom\ttstat\tpval_nominal\tslope\tslope_se\tpval_perm\tpval_beta"));
    thread thread1([&]() {
        startMPCset1(pheno_path, geno_path, pheno_pos, geno_pos, openPorts, 0, address, lo_row, mid_row, 0,permut,cislogger,nominallogger);
    });
    thread thread2([&]() {
        startMPCset1(pheno_path, geno_path, pheno_pos, geno_pos, openPorts, 12, address, mid_row, hi_row, 0, permut,cislogger2,nominallogger2);
    });

    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));

    for (auto& thread : threads) {
        thread.join();
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/60.0;
     if (getrusage(RUSAGE_SELF, &r_usage) == 0) {
        std::cout << "Memory usage: " << r_usage.ru_maxrss << " KB" << std::endl;
    } else {
        std::cerr << "Failed to get resource usage." << std::endl;
    }
    cout << "Execution time: " << totaldurationInminutes << " minutes" << endl;

    return 0;
}

// privateQTLII preprocessing
int main2(int argc, char* argv[])
{
    omp_set_num_threads(4);
    if (argc < 6)
    {
        cout << "Please provide at least: low, mid, high, zscore file, normalization method, and split set.\n";
        return 1;
    }
    int lo_row = stoi(argv[1]);
    int mid_row = stoi(argv[2]);
    int hi_row = stoi(argv[3]);
    string pheno_input = argv[4];
    string zscorefile = argv[5];
    string norm_method = argv[6];
    string output_file = argv[7];
    int siteA_n = stoi(argv[8]);
    int siteB_n = stoi(argv[9]);
    int siteC_n = stoi(argv[10]);

    int startingPort = 12300;
    int assignedPorts=0;
    vector<int> openPorts;
    while (assignedPorts < 24)
    {
        int port = startingPort;
        if(isPortOpen(port))
        {
            openPorts.push_back(port);
            assignedPorts++;
        }
        startingPort++;
    }
    vector<string> address{"localhost","localhost","localhost","localhost"};

    auto start = chrono::high_resolution_clock::now();
    int numCores = thread::hardware_concurrency();
    vector<thread> threads;
    vector<vector<double>> resultVec1, resultVec2;
    vector<string> gene_string1, gene_string2;
    thread thread1([&]() {
        startMPCset2(norm_method, pheno_input, openPorts, 0, address, lo_row, mid_row, zscorefile, 0,resultVec1,gene_string1);
    });
    thread thread2([&]() {
        startMPCset2(norm_method, pheno_input, openPorts, 12, address, mid_row, hi_row,zscorefile, 0,resultVec2,gene_string2);
    });

    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));

    for (auto& thread : threads) {
        thread.join();
    }
    resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    gene_string1.insert(gene_string1.end(), gene_string2.begin(), gene_string2.end());
    // writeNormalizedToTSV(resultVec1, gene_string1, "private_deseq2_invcdf_gmg_normalized");

    vector<vector<double>> siteA(resultVec1.size(), vector<double>(siteA_n));
    vector<vector<double>> siteB(resultVec1.size(), vector<double>(siteB_n));
    vector<vector<double>> siteC(resultVec1.size(), vector<double>(siteC_n));

    // Iterate through rows and columns to populate sets
    for (int i = 0; i < resultVec1.size(); ++i) {
        for (int j = 0; j < siteA_n; ++j) {
            siteA[i][j] = resultVec1[i][j];
        }

        for (int j = 0; j < siteB_n; ++j) {
            siteB[i][j] = resultVec1[i][siteA_n + j];
        }

        for (int j = 0; j < siteC_n; ++j) {
            siteC[i][j] = resultVec1[i][siteA_n + siteB_n + j];
        }
    }
    vector<vector<double>> siteA_cov, siteB_cov, siteC_cov;
    PCA(siteA, siteA_cov, 3);
    PCA(siteB, siteB_cov, 3);
    PCA(siteC, siteC_cov, 3);
    cout << "siteA_cov shape: " << siteA_cov.size() << ", " << siteA_cov[0].size() << endl;
    cout << "siteB_cov shape: " << siteB_cov.size() << ", " << siteB_cov[0].size() << endl;
    cout << "siteC_cov shape: " << siteC_cov.size() << ", " << siteC_cov[0].size() << endl;
    // writematrixToTSV(siteA_cov, "private_deseq2_invcdf_gmg_siteA_pc");
    // writematrixToTSV(siteB_cov, "private_deseq2_invcdf_gmg_siteB_pc");
    // writematrixToTSV(siteC_cov, "private_deseq2_invcdf_gmg_siteC_pc");
    Residualizer res1(siteA_cov);
    Residualizer res2(siteB_cov);
    Residualizer res3(siteC_cov);
    vector<vector<double>> res_A = res1.transform(siteA);
    vector<vector<double>> res_B = res2.transform(siteB);
    vector<vector<double>> res_C = res3.transform(siteC);
    for (int i = 0; i < res_A.size(); ++i) {
        res_A[i].insert(res_A[i].end(), res_B[i].begin(), res_B[i].end());
        res_A[i].insert(res_A[i].end(), res_C[i].begin(), res_C[i].end());
    } 
    cout << "residualized shape: \n" << res_A.size() << ", " << res_A[0].size() <<endl;
    // string to_write = "private_deseq2_invcdf_"+split_set+"_residualized_700"; //CHANGE
    writeNormalizedToTSV(res_A, gene_string1, output_file);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/60.0;
     if (getrusage(RUSAGE_SELF, &r_usage) == 0) {
        std::cout << "Memory usage: " << r_usage.ru_maxrss << " KB" << std::endl;
    } else {
        std::cerr << "Failed to get resource usage." << std::endl;
    }
    cout << "Total execution time: " << totaldurationInminutes << " minutes" << endl;

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Please provide at least one argument: privateQTL version.\n";
        return 1;
    }
    int newArgc = (argc - 1);
    char** newArgv = argv + 1;

    int version = stoi(argv[1]);
    if (version == 1) {
        cout << "=======Running privateQTL-I =======\n";
        return main1(newArgc, newArgv);
    } else if (version == 2) {
        cout << "=======Running privateQTL-II =======\n";
        return main2(newArgc, newArgv);
    } else {
        cout << "Invalid privateQTL version.\n";
        return 1;
    }
}