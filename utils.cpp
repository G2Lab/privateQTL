#include "utils.h"

vector<vector<double>> eigenMatrixToVectorOfVectors(const MatrixXd& input) {
    vector<vector<double>> result(input.rows(), vector<double>(input.cols()));
    for (int i = 0; i < input.rows(); ++i) {
        for (int j = 0; j < input.cols(); ++j) {
            result[i][j] = input(i, j);
        }
    }
    return result;
}

MatrixXd vectorOfVectorsToEigenMatrix(const vector<vector<double>>& input) {
    MatrixXd result; 
    result.resize(input.size(), input[0].size());
    // MatrixXd result(input.size(), input[0].size());
    for (size_t i = 0; i < input.size(); ++i) {
        for (size_t j = 0; j < input[i].size(); ++j) {
            result(i, j) = input[i][j];
        }
    }
    return result;
}
vector<double> CSVtoVector(const string &filename)
{
    vector<double> input_vec;
    // vector<string> string_vector;
    ifstream data(filename);
    string line;
    // int line_count = 0;
    // string cell;

    while (getline(data, line))
    {
        // string_vector.push_back(line);
        double entry = ::stof(line);
        input_vec.push_back(entry);
        // cout << weight << " ";
    }
    return input_vec;
}
vector<string> TSVtoVector(const string &filename) {
    vector<string> input_vec;
    ifstream data(filename);
    string line;

    while (getline(data, line)) {
        istringstream iss(line);
        string token;

        while (getline(iss, token, '\t')) {
            input_vec.push_back(token);
        }
    }

    return input_vec;
}
vector<double> TSVtoDoubleVector(const string &filename) {
    vector<double> input_vec;
    ifstream data(filename);
    string line;

    while (getline(data, line)) {
        istringstream iss(line);
        string token;

        while (getline(iss, token, '\t')) {
            double entry = stod(token);
            // rowVector.push_back(entry);
            input_vec.push_back(entry);
        }
    }
    return input_vec;
}
vector<double> getRowFromMatrixFile(string& filename, int rowIndex) {
    vector<double> rowVector;
    ifstream data(filename);
    string line;
    int currentRow = 0;
    while (getline(data, line)) {
        if (currentRow == 0) {
            // Skip the first row (header)
            currentRow++;
            continue;
        }
        if (currentRow-1 == rowIndex) {
            stringstream lineStream(line);
            string cell;
            int currentColumn = 0;

            while (getline(lineStream, cell, '\t')) 
            {
                if (currentColumn == 0) {
                    currentColumn++;
                    continue;
                }
                try {
                double entry = stod(cell);
                rowVector.push_back(entry);
                }
                catch (const exception& e) {
                    cerr << "Exception caught: " << e.what() << std::endl;
                }
                currentColumn++;
            }
        }
        currentRow++;
    }
    return rowVector;
}
void read_bedfile_row(vector<double>& rowData, string& geneID, const string& filename, int row, int skipcols, bool header) {
    // vector<double> rowData;
    ifstream data(filename);
    string line;
    int currentRow = 0;
    int index = 0;
    while (getline(data, line)) {
        if (header){
            if (currentRow == 0) {
            // Skip the first row (header)
            currentRow++;
            continue;
            }
            index = currentRow-1;
        }
        else {
            index = currentRow;
        }
        if (index == row) {
            stringstream lineStream(line);
            string cell;
            int currentColumn = 0;

            while (getline(lineStream, cell, '\t')) 
            {
                if (currentColumn == 0 && skipcols == 0) {
                    // cout << "currentColumn "<< currentColumn << endl;
                    geneID = cell;
                    currentColumn++;
                    continue;
                }
                else if (currentColumn < skipcols) {
                    if (currentColumn == skipcols-1)
                        geneID = cell;
                    currentColumn++;
                    continue;
                }
                try {
                double entry = stod(cell);
                rowData.push_back(entry);
                }
                catch (const exception& e) {
                    cerr << "Exception caught: " << e.what() << std::endl;
                }
                currentColumn++;
            }
        }
        currentRow++;
    }
    // // Skip the header if required
    // if (header)
    //     getline(data, line);
    // // Skip rows until the desired row
    // for (int currentRow = 0; currentRow < row-1; ++currentRow) {
    //     if (!getline(data, line)) {
    //         cerr << "Desired row not found." << endl;
    //     }
    // }
    // stringstream lineStream(line);
    // string cell;
    // // Skip the first N columns
    // for (int i = 0; i < skipcols; ++i) {
    //     if (!getline(lineStream, cell, '\t')) {
    //         cerr << "Not enough columns in the row." << endl;
    //     }
    //     if (i == skipcols-1)
    //         geneID = cell;

    // }
    // while (getline(lineStream, cell, '\t')) {
    //     try {
    //         double entry = stod(cell);
    //         rowData.push_back(entry);
    //     } catch (const exception& e) {
    //         cerr << "Exception caught: " << e.what() << endl;
    //     }
    // }
    data.close();
}
vector<vector<double>> getTPMFromMatrixFile(const string& filename, vector<string>& geneID) {
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        // int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, '\t');
        geneID.push_back(cell);
        getline(lineStream, cell, '\t');

        vector<double> rowVector;
        while (getline(lineStream, cell, '\t')) {
            try {
                double entry = stod(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            // currentColumn++;
        }

        rowsData.push_back(rowVector);
        // currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
vector<vector<double>> getCovariates(const string& filename) {
    // cout << "enetered getcovariates.\n";
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        // int currentColumn = 0;

        // Skip the first two columns
        getline(lineStream, cell, ',');
        // geneID.push_back(cell);
        // getline(lineStream, cell, '\t');

        vector<double> rowVector;
        while (getline(lineStream, cell, ',')) {
            try {
                double entry = stod(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            // currentColumn++;
        }

        rowsData.push_back(rowVector);
        // currentRow++;
    }

    // Close the file after reading
    data.close();
    cout << "leaving getcovariates.\n";
    return rowsData;
}
vector<vector<uint64_t>> getCountFromMatrixFile(const string& filename, vector<string>& geneID, int skipcols) {
    vector<vector<uint64_t>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        stringstream lineStream(line);
        string cell;
        // int currentColumn = 0;

        // Skip the first two columns
        for (int i=0; i<skipcols; i++){
            getline(lineStream, cell, '\t');
            if (i==0) {geneID.push_back(cell);}
        }
        // getline(lineStream, cell, '\t');
        // geneID.push_back(cell);
        // getline(lineStream, cell, '\t');

        vector<uint64_t> rowVector;
        while (getline(lineStream, cell, '\t')) {
            try {
                uint64_t entry = stoi(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cerr << "Exception caught: " << e.what() << endl;
            }
            // currentColumn++;
        }

        rowsData.push_back(rowVector);
        // currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}
vector<uint64_t> ScaleVector(vector<double> &v, int k)
{
    vector<uint64_t> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = v[i] * k;
    return intvec;
}
vector<vector<int64_t>> ScaleVector(vector<vector<double>>& v, int k)
{
    vector<vector<int64_t>> intvec(v.size(), vector<int64_t>(v[0].size()));
    for (int i=0; i<v.size(); i++)
    {
        for (int j=0; j<v[0].size(); j++)
        {
            intvec[i][j]=v[i][j]*k;
        }
    }
    return intvec;
}
vector<int64_t> ScaleVector_signed(vector<double> &v, int k)
{
    vector<int64_t> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = v[i] * k;
    return intvec;
}


vector<double> UnscaleVector_signed(vector<int64_t> &v, int k)
{
    vector<double> intvec(v.size(), 0);
    for (int i = 0; i < v.size(); ++i)
        intvec[i] = static_cast<double>(v[i]) / static_cast<double>(k);
    return intvec;
}
vector<uint64_t> ShiftVector(vector<int64_t>& vec, uint64_t number) 
{
    vector<uint64_t> shifted(vec.size());
    int64_t shiftedValue;
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i]<0)
            shiftedValue = vec[i] + number;
        else
            shiftedValue = vec[i];
        shifted[i] = static_cast<uint64_t>(shiftedValue);
    }
    return shifted;
}

vector<vector<uint64_t>> ShiftVector(vector<vector<int64_t>>& vec, uint64_t number) 
{
    vector<vector<uint64_t>> shifted(vec.size(), vector<uint64_t>(vec[0].size()));
    int64_t shiftedValue;
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j=0; j< vec[0].size();j++)
        {
            if(vec[i][j]<0)
            {
                shiftedValue = vec[i][j] + number;
            }
            else{
                shiftedValue = vec[i][j];
            }
            shifted[i][j] = static_cast<uint64_t>(shiftedValue);
        }
        
    }
    return shifted;
}

vector<int64_t> UnshiftVector(vector<uint64_t>& shifted, uint64_t number)
{
    vector<int64_t> unshifted(shifted.size());
    for (size_t i = 0; i < shifted.size(); ++i) {
        if (shifted[i] > number/2)
            unshifted[i] = static_cast<int64_t>(shifted[i]) - number;
        else
            unshifted[i] = static_cast<int64_t>(shifted[i]);
        // unshifted[i] = static_cast<int64_t>(shiftedValue);
    }
    return unshifted;
}

template <typename T>
void print_vector(vector<vector<T>> printme)
{
    for (int i =0; i<printme.size(); i++)
    {
        print_vector(printme[i]);
    }
    std::cout << "\n";
}

template <typename T>
void print_vector(Vec<T> printme)
{
    for (int i =0; i<printme.length(); i++)
    {
        std::cout << printme.get(i) << " ";
    }
    std::cout << "\n";
}
uint64_t nearestPowerOf2(int N)
{
    uint64_t a = log2(N);
 
    // if (pow(2, a) == N)
    //     return N;
 
    return pow(2, a + 1);
}
double betaLogLikelihood(const gsl_vector *v, void *params) {
	double * p = (double *) params;
	double beta_shape1 = gsl_vector_get(v, 0);
	double beta_shape2 = gsl_vector_get(v, 1);
	return -1.0 * ((beta_shape1 - 1) * p[0] + (beta_shape2 - 1) * p[1] - p[2] * gsl_sf_lnbeta(beta_shape1, beta_shape2));
}

int mleBeta(vector < double > & pval, double & beta_shape1, double & beta_shape2) {

	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, beta_shape1);
	gsl_vector_set (x, 1, beta_shape2);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (2);
	gsl_vector_set (ss, 0, beta_shape1/10);
	gsl_vector_set (ss, 1, beta_shape2/10);

	//Initialize method and iterate
	double par [3];
	par[0] = 0.0;
	par[1] = 0.0;
	for (int e = 0 ; e < pval.size(); e ++) {
		if (pval[e] == 1.0) pval[e] = 0.99999999;
		par[0] += log (pval[e]);
		par[1] += log (1 - pval[e]);
	}
	par[2] = pval.size();
	gsl_multimin_function minex_func;
	minex_func.n = 2;
	minex_func.f = betaLogLikelihood;
	minex_func.params = par;

	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	//Optimization iteration
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
	} while (status == GSL_CONTINUE && iter < 1000);

	//Output new beta shape values
	beta_shape1 = gsl_vector_get (s->x, 0);
	beta_shape2 = gsl_vector_get (s->x, 1);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}

double doublemean(vector < double > & X) {
		double mean = 0.0;
		for (int x = 0 ; x < X.size() ; x ++) mean += X[x];
		mean /= X.size();
		return mean;
	}

double doublevariance(vector < double > & X, double mean) {
    double variance = 0.0;
    for (int x = 0 ; x < X.size() ; x++) variance += (X[x] - mean) * (X[x] - mean);
    variance /= (X.size() - 1);
    return variance;
}
struct data_to_function {
	int n;
	double * C;

	data_to_function (int _n, double * _C) {
		n = _n;
		C = _C;
	}

};
double degreeOfFreedom(const gsl_vector *v, void *params) {
    // pair<vector<double>*, double>* p = static_cast<pair<vector<double>*, double>*>(params);	
    // vector<double>& corr = *(p->first);
    // vector < double > pval = vector < double >(corr.size(), 0.0);
    // pair<int, double*>* p = static_cast<pair<int, double*>*> (params);
    data_to_function * d = (data_to_function *) params;
    // int size = p->first;
    // double* corrs = p->second;
    vector < double > pval = vector < double >(d->n, 0.0);
	double mean = 0.0;
	for (int c = 0 ; c < d->n ; c++) {
		pval[c] = getPvalue(d->C[c], gsl_vector_get(v, 0));
		mean += pval[c];
	}
	mean /= pval.size();
	double variance = 0.0;
	for (int c = 0 ; c < pval.size() ; c++) variance += (pval[c] - mean) * (pval[c] - mean);
	variance /= (pval.size() - 1);

	double shape2 = abs((mean * (mean * (1 - mean ) / variance - 1)) - 1.0);
	//cout << "O = " << mean << " " << shape2 << endl;
	return shape2;
}
int learnDF(vector < double > & corr, double & df) {
    // writeVectorToTSV(corr,string("testr_perm"));
	//Set starting point to moment matching estimates
	gsl_vector * x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, df);

	//Set initial step sizes to shape1 and shape2 scales
	gsl_vector * ss = gsl_vector_alloc (1);
	gsl_vector_set (ss, 0, df * 0.1);

	// pair<vector<double>*, double> params(&corr, df);
    data_to_function * par  = new data_to_function (corr.size(), &corr[0]);
    // cout << "1" << endl;
    // pair<int, double*> params (corr.size(), &corr[0]);
    // cout << "2" << endl;
	gsl_multimin_function minex_func;
	minex_func.n = 1;
	minex_func.f = degreeOfFreedom;
	// minex_func.params = &params;
    minex_func.params = (void*)par;
    // cout << "3" << endl;
	//Initialize optimization machinery
	const gsl_multimin_fminimizer_type * T = gsl_multimin_fminimizer_nmsimplex2;
    // cout << "3.1" << endl;
	gsl_multimin_fminimizer * s = gsl_multimin_fminimizer_alloc (T, 1);
    // cout << "3.2" << endl;
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    // cout << "4" << endl;
	//Optimization iteration
	//cout << "\n ========================" << endl;
	size_t iter = 0;
	int status;
	double size;
	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status) break;
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 0.01);
		// printf ("%d %10.3e f() = %7.10f size = %.10f\n", iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < 20);

	//Output new beta shape values
	df = gsl_vector_get (s->x, 0);

	//Free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	return (status == GSL_SUCCESS);
}
double getTstat2(double corr, double df) {
    return df * corr * corr / (1 - corr * corr);
}
double getPvalueFromTstat2(double tstat2, double df) {
    return pf(tstat2, 1, df, 0, 0);
}
double getPvalue(double corr, double df) {
    assert(corr != 1 && corr != -1);
    assert(df * corr * corr != (1 - corr * corr));
    return pf(df * corr * corr / (1 - corr * corr), 1, df,0,0);
}
double getSlope(double nominal_correlation, double psd, double gsd) {
    if (gsd < 1e-16 || psd < 1e-16) return 0;
    else return nominal_correlation * psd / gsd;
}

void svd_flip(MatrixXd& u, MatrixXd& v, bool u_based_decision) {
    std::cout << "entered svd_flip." << std::endl;
    if (u_based_decision) {
        // Columns of u, rows of v
        VectorXi max_abs_cols(u.cols());

        for (int j = 0; j < u.cols(); ++j) {
            double max_abs_value = 0.0;
            int max_abs_index = 0;

            for (int i = 0; i < u.rows(); ++i) {
                double abs_value = std::abs(u(i, j));
                if (abs_value > max_abs_value) {
                    max_abs_value = abs_value;
                    max_abs_index = i;
                }
            }

            max_abs_cols(j) = max_abs_index;
        }
        // recreating square matrix with max_abs_cols
        MatrixXd u_square(u.cols(), u.cols());
        for (int i = 0; i < u.cols(); ++i) {
            u_square.row(i) = u.row(max_abs_cols(i));
        }
        // Creating signs vector based on signs of diagonal element 
        VectorXd signs(u.cols());
        for (int i = 0; i < u.cols(); ++i) {
            signs(i) = (u_square.diagonal()(i) >= 0) ? 1.0 : -1.0;
        }
        
        for (int i = 0; i < u.cols(); ++i) {
            u.col(i) *= signs(i);
        }
        // u.array() *= signs.array();
        for (int i=0; i<v.rows(); ++i){
            v.row(i) *= signs(i);
        }
        // v.array().colwise() *= signs.col(0).array();
    } else {
        // Rows of v, columns of u
        VectorXi max_abs_rows = (v.array().abs()).rowwise().maxCoeff().cast<int>();

        MatrixXd signs(v.rows(), v.cols());
        for (int i = 0; i < v.rows(); ++i) {
            signs.row(i) = (v.row(i).array() >= 0).template cast<double>() * 2.0 - 1.0;
        }

        u.array().rowwise() *= signs.row(0).array();
        v.array() *= signs.array();
    }
}
void PCA(vector<vector<double>>& data, vector<vector<double>>& pc, int n_components)
{
    MatrixXd data_mat = vectorOfVectorsToEigenMatrix(data);
    VectorXd mean = data_mat.colwise().mean();
    data_mat.rowwise() -= mean.transpose();
    // Compute covariance matrix
    MatrixXd covMatrix = (data_mat.transpose() * data_mat) / (data_mat.rows() - 1);
    if (data.size() < 500 && data[0].size() < 500){
        cout << "Eigen.." << endl;
        // Perform eigendecomposition
        Eigen::EigenSolver<MatrixXd> solver(covMatrix);
        VectorXd eigenvalues = solver.eigenvalues().real();

        // VectorXd explained_variance_ratio = eigenvalues / eigenvalues.sum();
        // VectorXd cumulative_variance(explained_variance_ratio.size());
        // double sum = 0.0;
        // for (int i = 0; i < explained_variance_ratio.size(); ++i) {
        //     sum += explained_variance_ratio(i);
        //     cumulative_variance(i) = sum;
        // }
        // // Find the optimal number of components
        // int optimal_components = 0;
        // while (optimal_components < cumulative_variance.size() &&
        //        cumulative_variance[optimal_components] < cumulative_variance_threshold) {
        //     optimal_components++;
        // }
        // cout << "Optimal components: " << optimal_components << endl;
        MatrixXd eigenvectors = solver.eigenvectors().real();
        eigenvectors = eigenvectors.leftCols(n_components);
        pc = eigenMatrixToVectorOfVectors(eigenvectors);
    }
    else {
        cout << "Randomized" << endl;
        Eigen::JacobiSVD<MatrixXd> svd(covMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
        VectorXd singularValues = svd.singularValues();
        // reverse(singularValues.data(), singularValues.data() + singularValues.size()); // Reverse the singular values
        // cout << singularValues(0) << ", " <<singularValues(1) << ", " << singularValues(2) << endl;
        // VectorXd explained_variance_svd = singularValues.array().square() / (singularValues.size() - 1);
        // double total_var_svd = explained_variance_svd.sum();
        // VectorXd explained_variance_ratio_svd = explained_variance_svd / total_var_svd;
        // // VectorXd explained_variance_ratio_svd = singularValues.array().square() / (singularValues.array().square().sum());

        // VectorXd cumulative_variance_svd(explained_variance_ratio_svd.size());
        // double sum = 0.0;
        // for (int i = 0; i < explained_variance_ratio_svd.size(); ++i) {
        //     sum += explained_variance_ratio_svd(i);
        //     cumulative_variance_svd(i) = sum;
        // }
        // // Find the optimal number of components for SVD
        // int optimal_components_svd = 0;
        // while (optimal_components_svd < cumulative_variance_svd.size() &&
        //        cumulative_variance_svd[optimal_components_svd] < cumulative_variance_threshold) {
        //     optimal_components_svd++;
        // }
        // if (optimal_components_svd == 0){
        //     optimal_components_svd = 1;
        // }
        MatrixXd u_copy = svd.matrixU();
        MatrixXd v_copy = svd.matrixV();
        svd_flip(u_copy, v_copy, true);
        MatrixXd components = v_copy.leftCols(n_components);
        pc = eigenMatrixToVectorOfVectors(components);
    }
}
vector<double> center_normalize(vector<vector<double>>& M) {
    int rows = M.size();
    int cols = M[0].size();

    // Compute the mean and variance for each row
    vector<vector<double>> N(rows, vector<double>(cols, 0.0));
    vector<double> row_variances(rows, 0.0);
    for (int i = 0; i < rows; ++i) {
        // Calculate row mean
        double row_mean = 0.0;
        for (int j = 0; j < cols; ++j) {
            row_mean += M[i][j];
        }
        row_mean /= cols;

        // Center and normalize the row, and calculate variance
        for (int j = 0; j < cols; ++j) {
            N[i][j] = M[i][j] - row_mean;
            row_variances[i] += N[i][j] * N[i][j];
        }
    }

    // Normalize each row by its variance
    for (int i = 0; i < rows; ++i) {
        double norm = sqrt(row_variances[i]);
        for (int j = 0; j < cols; ++j) {
            N[i][j] /= norm;
        }
        row_variances[i] /= (cols-1); // Not the same
    }
    swap(N, M);
    return row_variances;
}
double center_normalize_vec(vector<double>& row) {
    int cols = row.size();

    // Calculate the mean of the row
    double row_mean = 0.0;
    for (int j = 0; j < cols; ++j) {
        row_mean += row[j];
    }
    row_mean /= cols;

    // Center and normalize the row, and calculate variance
    double row_variance = 0.0;
    for (int j = 0; j < cols; ++j) {
        row[j] -= row_mean;
        row_variance += row[j] * row[j];
    }

    // Normalize the row by its variance
    double norm = sqrt(row_variance);
    for (int j = 0; j < cols; ++j) {
        row[j] /= norm;
    }
    row_variance /= (cols - 1);

    return row_variance;
}

// template <typename T>
// void writematrixToTSV(const vector<vector<T>>& data, const string& name)
// {
//     string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name  + ".tsv";
//     ofstream file(filename);
//     if (file.is_open())
//     {
//         for (int i = 0; i < data.size(); ++i)
//         {
//             for (const T& value : data[i])
//             {
//                 file << value << "\t";
//             }
//             file << std::endl;
//         }
//         file.close();
//         cout << string(name+" matrix successfully written to TSV file.") << endl;
//     }
//     else
//     {
//         cout << "Error opening the file." << endl;
//     }
// } 
void writeNormalizedToTSV(const vector<vector<double>>& data, const vector<string>& gene_strings, const string& name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name  + ".tsv";
    ofstream file(filename);
    if (file.is_open())
    {
        for (int i = 0; i < data.size(); ++i)
        {
            // Add gene name in the first column
            file << gene_strings[i] << "\t";

            // Add the rest of the values in the row
            for (const double& value : data[i])
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


vector<vector<double>> getMatrixFile(const string& filename, int startrow, int endrow, bool header, bool index) {
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    int currentRow = 0;

    if (header)// Skip the first row (header)
        getline(data, line);

    while (getline(data, line) && currentRow < endrow) {
        if (currentRow >= startrow) { // Start reading from startrow
            stringstream lineStream(line);
            string cell;
            int currentColumn = 0;

            // Skip the first column
            if (index)
                getline(lineStream, cell, '\t');

            vector<double> rowVector;
            while (getline(lineStream, cell, '\t')) {
                try {
                    double entry = stod(cell);
                    rowVector.push_back(entry);
                } catch (const exception& e) {
                    cerr << "Exception caught: " << e.what() << endl;
                }
                currentColumn++;
            }

            rowsData.push_back(rowVector);
        }

        currentRow++;
    }

    // Close the file after reading
    data.close();

    return rowsData;
}

