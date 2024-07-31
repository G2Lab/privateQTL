#ifndef PRIVATEQTL_CP_H
#define PRIVATEQTL_CP_H
// #define EIGEN_USE_THREADS
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/Defines.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <map>
#include <set>
#include <bitset>
#include <fstream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm> //for debugging
#include <sys/socket.h>
#include <netinet/in.h>
#include <cstring>
#include <Eigen/Dense>
#include "utils.h"
#include <omp.h>
// using namespace Eigen;
using namespace osuCrypto;
using namespace std;
using namespace NTL;
// namespace R = RmathNamespace;
class Logger {
public:
    Logger(const std::string& filename) {
        logFile.open(filename, std::ios::out | std::ios::trunc);
    }

    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    void log(const std::string& message) {
        if (logFile.is_open()) {
            std::time_t currentTime = std::time(nullptr);
            logFile << message << std::endl;
        }
    }

private:
    std::ofstream logFile;
};
class mpc
{
public:
    IOService ios;
    PRNG globalprng;
    int pid;
    uint64_t inv;
    uint64_t lk;
    uint64_t n;
    uint64_t p;
    double shiftsize;
    typedef Eigen::Matrix <uint64_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXi;
    
    mpc(){};
    bool initialize(int pid, string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2);
    bool setupSeeds(); // shared seeds for PRNG
    bool areVectorsEmpty() const {
        return shares.empty();// && identity.empty() && zscores.empty();
    }
    void ready();
    void receiveSecrets();
    void assertSize(vector<ZZ_p> pi, string tag);
    void assertValid(vector<ZZ_p> sigma, string tag);
    vector<uint64_t> apply_plaintext_perm(vector<uint64_t> rho, vector<uint64_t> sigma);
    vector<ZZ_p> Frand(uint64_t bufferSize);
    vector<ZZ_p> reveal(vector<ZZ_p> &pi, bool isperm);
    void reshare(vector<ZZ_p> &shares, int reshareID);
    void reshareM(vector<vector<ZZ_p>>& shares, int reshareID);
    vector<ZZ_p> Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i);
    vector<ZZ_p> genbitperm(vector<ZZ_p> &keybit);
    vector<ZZ_p> inversePerm(vector<ZZ_p> pi);
    void apply_perm_local(bool participate,vector<ZZ_p> &v, vector<ZZ_p> &pi);
    void apply_perm_localM(bool participate, vector<vector<ZZ_p>> &v, vector<ZZ_p> &pi);
    void center_normalize_pheno();
    void permutPheno(int permut);
    void shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a);
    void unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b);
    void reveal_matrix(vector<vector<ZZ_p>>& geno, vector<vector<ZZ_p>>& pheno,string name);
    void calc_corr(Logger& cislogger, Logger& nominalLogger);
    void shuffleM(vector<ZZ_p> &pi, vector<vector<ZZ_p>> &a);
    void apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k);
    void compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho);
    // vector<ZZ_p> get_shared_inverse(vector<ZZ_p> sigma);
    void genperm(int row, string norm_method);
    void clearVectors() {
        shares.clear();
        for (auto &innerVec : geno) {
        innerVec.clear(); 
        }
        for (auto &innerVec : permutMat) {
        innerVec.clear(); 
        }
        geno.clear();
        permutMat.clear();
    }
    void receivePheno();
    void logRatio();
    void center_normalize_geno();
    void receiveGeno();
    vector<vector<ZZ_p>> matmult(vector<vector<ZZ_p>>& mat1, vector<vector<ZZ_p>>& mat2);
    void close();
    int validgene();
    inline void ZZtoEigen(vector<vector<ZZ_p>>& v, MatrixXi& dest1, MatrixXi& dest2, bool transposed) {
        int numRowsV = v.size();
        int numColsV = v[0].size();

        if (transposed) {
            // Check dimensions
            if (dest1.rows() != numColsV / 2 || dest1.cols() != numRowsV ||
                dest2.rows() != numColsV / 2 || dest2.cols() != numRowsV) {
                throw invalid_argument("Set eigen matrix size to be equal, please.");
            }

            for (int i = 0; i < numRowsV; i++) {
                for (int j = 0; j < numColsV / 2; j++) {
                    dest1(j, i) = conv<uint64_t>(v[i][2 * j]);
                    dest2(j, i) = conv<uint64_t>(v[i][2 * j + 1]);
                }
            }
        } else {
            // Check dimensions
            if (dest1.rows() != numRowsV || dest1.cols() != numColsV / 2 ||
                dest2.rows() != numRowsV || dest2.cols() != numColsV / 2) {
                throw invalid_argument("Set eigen matrix size to be equal, please.");
            }

            for (int i = 0; i < numRowsV; i++) {
                for (int j = 0; j < numColsV / 2; j++) {
                    dest1(i, j) = conv<uint64_t>(v[i][2 * j]);
                    dest2(i, j) = conv<uint64_t>(v[i][2 * j + 1]);
                }
            }
        }
    }


    inline void EigentoZZ(vector<uint64_t>& share1, vector<uint64_t>& share2, vector<vector<ZZ_p>>& dest){
        // vector<uint32_t> converted(v.size());
        if (share1.size() != share2.size())
            throw invalid_argument("Your shares have different sizes.");
        if (share1.size() != dest.size()*dest[0].size()/2)
            throw invalid_argument("Your shares and destination matrix have different sizes.");
        for (int i=0; i<dest.size(); i++)
        {
            for (int j=0; j<dest[0].size()/2; j++)
            {
                dest[i][2*j] = conv<ZZ_p>(share1[i*dest[0].size()/2+j]);
                dest[i][2*j+1] = conv<ZZ_p>(share2[i*dest[0].size()/2+j]);
            }
        }
    }
    void writeEigenToTSV(Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat, const string& name);
private:
    block commonSeed = oc::toBlock(27);
    map<int, PRNG *> seedpair;
    PRNG* localprng;
    vector<ZZ_p> shares;
    vector<vector<ZZ_p>> permutMat;
    vector<ZZ_p> identity;
    vector<ZZ_p> zscores;
    vector<vector<ZZ_p>> geno;
    vector<vector<ZZ_p>> pheno;
    vector<uint64_t> shape;
    // atomic<int>& readyCounter;
    // atomic<bool> allPartiesReady;
    Channel dataowner;
    Channel toOwner;
    Channel fromPlus;
    Channel fromMinus;
    Channel toPlus;
    Channel toMinus;  
};
inline vector<ZZ_p> convVec(vector<uint64_t> v){
    vector<ZZ_p> converted(v.size());
    for (int i=0; i<v.size(); i++)
    {
        converted[i] = to_ZZ_p(v[i]);
    }
    return converted;
}
inline vector<uint64_t> convVec(vector<ZZ_p> v){
    vector<uint64_t> converted(v.size());
    for (int i=0; i<v.size(); i++)
    {
        converted[i] = conv<uint64_t>(v[i]);
    }
    return converted;
}
template <typename T>
void print_vector(vector<T> printme)
{
    for (int j = 0; j < printme.size(); j++)
    {
        std::cout << printme[j] << " ";
    }
    std::cout << "\n";
}



template <typename T>
void assert_equal(vector<T> one, vector<T> two, string tag)
{
    if (one.size() != two.size()) {
        throw logic_error("[" + tag + "]: vectors are not equal length");
    }
    for (int i = 0; i < one.size(); i++)
    {
        if (one[i] != two[i]) {
            throw logic_error("[" + tag + "]: vectors are not equal");
        }
    }
}
// inline double variance(vector < double > & X, double mean) {
//     double variance = 0.0;
//     for (int x = 0 ; x < X.size() ; x++) variance += (X[x] - mean) * (X[x] - mean);
//     variance /= (X.size() - 1);
//     return variance;
// }

// inline double mean(vector < double > & X) {
//     double mean = 0.0;
//     for (int x = 0 ; x < X.size() ; x ++) mean += X[x];
//     mean /= X.size();
//     return mean;
// }


// void ZZtoEigen(vector<vector<ZZ_p>>& v, MatrixXd eigenMatrix& dest){
//     // vector<uint32_t> converted(v.size());
//     if (dest.rows() != v.size() | dest.cols() != v[0].size())
//         throw invalid_argument("Set eigen matrix size to be equal, please.");
//     for (int i=0; i<v.size(); i++)
//     {
//         for (int j=0; j<v[0].size(); j++)
//         {
//             dest(i,j) = conv<uint32_t>(v[i][j]);
//         }
//     }
// }

// void EigentoZZ(MatrixXd eigenMatrix& v, vector<vector<ZZ_p>>& dest){
//     // vector<uint32_t> converted(v.size());
//     if (v.rows() != dest.size() | v.cols() != dest[0].size())
//         throw invalid_argument("Set ZZ_p matrix size to be equal, please.");
//     for (int i=0; i<v.rows(); i++)
//     {
//         for (int j=0; j<v.cols(); j++)
//         {
//             dest[i][j] = conv<ZZ_p>(v(i,j));
//         }
//     }
// }

#endif