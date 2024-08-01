#ifndef UTILS_H
#define UTILS_H
#define MATHLIB_STANDALONE
#define EIGEN_USE_THREADS
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <Rmath.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <NTL/vec_ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#define BETA_SHAPE1_MIN 0.1
#define BETA_SHAPE2_MIN 1
#define BETA_SHAPE1_MAX 10
#define BETA_SHAPE2_MAX 1000000	
using namespace std;
using namespace Eigen;
using namespace NTL;
typedef Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixDR;
// Helper function to convert vector of vectors to Eigen matrix
// MatrixXd vectorOfVectorsToEigenMatrix(const vector<vector<double>>& input) {
//     MatrixXd result; 
//     result.resize(input.size(), input[0].size());
//     // MatrixXd result(input.size(), input[0].size());
//     for (size_t i = 0; i < input.size(); ++i) {
//         for (size_t j = 0; j < input[i].size(); ++j) {
//             result(i, j) = input[i][j];
//         }
//     }
//     return result;
// }
// Helper function to convert Eigen matrix to vector of vectors
// vector<vector<double>> eigenMatrixToVectorOfVectors(const MatrixXd& input) {
//     vector<vector<double>> result(input.rows(), vector<double>(input.cols()));
//     for (int i = 0; i < input.rows(); ++i) {
//         for (int j = 0; j < input.cols(); ++j) {
//             result[i][j] = input(i, j);
//         }
//     }
//     return result;
// }
vector<vector<double>> eigenMatrixToVectorOfVectors(const MatrixXd& input); 
MatrixXd vectorOfVectorsToEigenMatrix(const vector<vector<double>>& input);
vector<double> CSVtoVector(const string &filename);
vector<string> TSVtoVector(const string &filename);
vector<double> TSVtoDoubleVector(const string &filename);
vector<double> getRowFromMatrixFile(string& filename, int rowIndex);
void read_bedfile_row(vector<double>& rowData, string& geneID, const string& filename, int row, int skipcols, bool header);
vector<vector<double>> getTPMFromMatrixFile(const string& filename, vector<string>& geneID);
vector<vector<double>> getCovariates(const string& filename);
vector<vector<uint64_t>> getCountFromMatrixFile(const string& filename, vector<string>& geneID, int skipcols);
vector<uint64_t> ScaleVector(vector<double> &v, int k);
vector<vector<int64_t>> ScaleVector(vector<vector<double>>& v, int k);
vector<int64_t> ScaleVector_signed(vector<double> &v, int k);
vector<double> UnscaleVector_signed(vector<int64_t> &v, int k);
vector<uint64_t> ShiftVector(vector<int64_t>& vec, uint64_t number); 
vector<vector<uint64_t>> ShiftVector(vector<vector<int64_t>>& vec, uint64_t number);
vector<int64_t> UnshiftVector(vector<uint64_t>& shifted, uint64_t number);
uint64_t nearestPowerOf2(int N);
int mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2);
double doublemean(vector < double > & X);
double doublevariance(vector < double > & X, double mean);
double degreeOfFreedom(const gsl_vector *v, void *params);
int learnDF(vector < double > & corr, double & df);
double getTstat2(double corr, double df);
double getPvalueFromTstat2(double tstat2, double df);
double getPvalue(double corr, double df);
double getSlope(double nominal_correlation, double psd, double gsd);
void svd_flip(MatrixXd& u, MatrixXd& v, bool u_based_decision);
void PCA(vector<vector<double>>& data, vector<vector<double>>& pc, int n_components);
vector<double> center_normalize(vector<vector<double>>& M);
double center_normalize_vec(vector<double>& M);
// template <typename T>
// void writematrixToTSV(const vector<vector<T>>& data, const string& name);
template <typename T>
void writematrixToTSV(const vector<vector<T>>& data, const string& name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name  + ".tsv";
    ofstream file(filename);
    if (file.is_open())
    {
        for (int i = 0; i < data.size(); ++i)
        {
            for (const T& value : data[i])
            {
                file << value << "\t";
            }
            file << std::endl;
        }
        file.close();
        cout << string(name+" matrix successfully written to TSV file.") << endl;
    }
    else
    {
        cout << "Error opening the file." << endl;
    }
}
void writeNormalizedToTSV(const vector<vector<double>>& data, const vector<string>& gene_strings, const string& name);
template <typename T>
void writeVectorToTSV(const vector<T>& data, string name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + ".tsv";
    ofstream file(filename);
    if (file.is_open())
    {
        for (const T& value : data)
        {
            file << value << "\t";
        }
        file.close();
        cout << string(name+ " vector successfully written to CSV file.") << endl;
    }
    else
    {
        cout << "Error opening the file." << endl;
    }
} 
vector<vector<double>> getMatrixFile(const string& filename, int startrow, int endrow, bool header, bool index);
template <typename T>
void print_vector(vector<vector<T>> printme);
template <typename T>
void print_vector(Vec<T> printme);
class Residualizer {
public: 
    Residualizer(const vector<vector<double>>& covariates) {
        // Convert vector of vectors to Eigen matrix
        MatrixXd C_t = vectorOfVectorsToEigenMatrix(covariates);
        cout << "inside constructor: "<< covariates.size() << ", " << covariates[0].size() << endl;
        Eigen::VectorXd col_means = C_t.colwise().mean();
        cout << col_means << endl;
        MatrixXd centered_C_t = C_t.rowwise() - col_means.transpose();
        Eigen::HouseholderQR<MatrixXd> qr(centered_C_t);
        // Q_t = qr.householderQ();
        Q_t = qr.householderQ() * MatrixXd::Identity(centered_C_t.rows(),centered_C_t.cols());
        dof = C_t.rows() - 2 - C_t.cols();
        // cout << centered_C_t(0,0) << ", " << centered_C_t(0,1) << ", " << centered_C_t(0,2) << endl;
        // cout << centered_C_t(1,0) << ", " << centered_C_t(1,1) << ", " << centered_C_t(1,2) << endl;
    }
    vector<vector<double>> transform(const vector<vector<double>>& M_t, bool center = true) {
        // Convert vector of vectors to Eigen matrix
        cout << "start transform function. "<< endl;
        MatrixXd M_t_matrix = vectorOfVectorsToEigenMatrix(M_t);
        
        // Residualize rows of M wrt columns of C
        cout << "befpre M0_t" << endl;
        MatrixXd M0_t = M_t_matrix.colwise() - M_t_matrix.rowwise().mean();
        cout << "M0_t made." << endl;
        if (center) {
            M0_t = M0_t - M0_t * Q_t * Q_t.transpose();
        } else {
            M0_t = M_t_matrix - M0_t * Q_t * Q_t.transpose();
        }

        // Convert Eigen matrix to vector of vectors
        return eigenMatrixToVectorOfVectors(M0_t);
    }
    vector<double> transform(const vector<double>& M_t, bool center = true) {
    // Convert vector to Eigen row matrix
        // cout << "transform entered." << endl;
        assert(!M_t.empty()); 
        Eigen::Map<const Eigen::Matrix<double, 1, Eigen::Dynamic>> M_t_matrix(M_t.data(), 1, M_t.size());

        // Residualize row of M wrt columns of C
        Eigen::Matrix<double, 1, Eigen::Dynamic> M0_t = M_t_matrix.array() - M_t_matrix.mean();
        // cout << M0_t.rows()<<","<<M0_t.cols() << "/"<< Q_t.rows() << "," << Q_t.cols() << endl;
        if (center) {
            M0_t.noalias() = M0_t - (M0_t * Q_t) * Q_t.transpose();
        } else {
            M0_t.noalias() = M_t_matrix - (M0_t * Q_t) * Q_t.transpose();
        }

        // Convert Eigen matrix to vector
        std::vector<double> result(M0_t.data(), M0_t.data() + M0_t.cols());
        return result;
    }
private:
    MatrixXd Q_t;
    int dof;
};
#endif
