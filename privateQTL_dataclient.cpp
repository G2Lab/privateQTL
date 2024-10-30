
#include "privateQTL_cp.h"
#include "privateQTL_dataclient.h"
#include <cmath>
#include "utils.h"
#include <chrono>

#include <sys/resource.h>
void dataclient::initialize(int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3) 
{
    // IOService ios;
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Endpoint epsend3(ios, address3, sendport3, EpMode::Server);
    Endpoint eprecv1(ios, address1, recvport1, EpMode::Client);
    Endpoint eprecv2(ios, address2, recvport2, EpMode::Client);
    Endpoint eprecv3(ios, address3, recvport3, EpMode::Client);
    this->owner_p1 = epsend1.addChannel();
    this->owner_p2 = epsend2.addChannel();
    this->owner_p3 = epsend3.addChannel();
    this->p1_owner = eprecv1.addChannel();
    this->p2_owner = eprecv2.addChannel();
    this->p3_owner = eprecv3.addChannel();
    
    NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
    p = pow(2,50);
    ZZ_p::init(to_ZZ(p)); 
    // cout << "P: " << p << endl;
    owner_p1.send(p);
    owner_p2.send(p);
    owner_p3.send(p);
    cout << "Owner established channels with the computing parties.\n";
}

prepareInput dataclient::load_genotype(string pheno_pos, string geno_path, string geno_pos, int cis_window)
{
    cout << "Loading Genotype matrix..." << flush;
    auto loadgeno = chrono::high_resolution_clock::now();
    prepareInput testinput(pheno_pos, geno_path ,geno_pos,cis_window);
    auto loadend = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = loadend - loadgeno;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    cout << durationInminutes << " minutes" << endl;
    return testinput;   
}
void dataclient::do_recv_string(vector<string>& stringvec)
{
    // this->dataowner.recv(snps);
    int buffersize;
    this->p1_owner.recv(buffersize);
    vector<char> serialized(buffersize);
    this->p1_owner.recv(serialized.data(), buffersize);
    string serializedVar(serialized.data(), buffersize);
    // vector<string> cisVar;
    istringstream iss(serializedVar);
    string token;
    while (getline(iss, token, ';')) {
        stringvec.push_back(token);
    }
}
void dataclient::do_send_string(vector<string>& tosend)
{
    string serializedvariants="";
    for (const std::string& str : tosend) {
        serializedvariants += str + ";"; // Use a suitable delimiter
    }
    int estimatedSize = static_cast<int>(serializedvariants.size());
    int bufferSize = estimatedSize;
    // Send the serialized data over the channel
    this->owner_p1.send(bufferSize);
    this->owner_p1.send(serializedvariants.data(), serializedvariants.size());
}
void dataclient::find_common(string phenofile)
{
    vector<vector<double>> common_pheno, common_geno;
    vector<string> snpset, geneset;
    do_send_string(this->myinput.snpIDs);
    // vector<string> snpset;
    do_recv_string(snpset);
    // p1_owner.recv(snpset);
    vector<string> geneIDs;
    vector<vector<double>> phenotype = getTPMFromMatrixFile(phenofile,geneIDs, false);
    do_send_string(geneIDs);
    // owner_p1.send(geneIDs);
    // vector<string> geneset;
    do_recv_string(geneset);
    // p1_owner.recv(geneset);
    unordered_set<string> geneset_set(geneset.begin(), geneset.end());
    unordered_set<string> snpset_set(snpset.begin(), snpset.end());
    for (size_t i = 0; i < geneIDs.size(); ++i) {
        if (geneset_set.find(geneIDs[i]) != geneset_set.end()) {
            common_pheno.push_back(phenotype[i]);
        }
    }
    cout << "pheno size: " << common_pheno.size() <<","<<common_pheno[0].size()<< endl;
    swap(this->pheno, common_pheno);
    swap(this->genes, geneset);
    for (size_t i = 0; i < myinput.snpIDs.size(); ++i) {
        if (snpset_set.find(myinput.snpIDs[i]) != snpset_set.end()) {
            common_geno.push_back(myinput.geno[i]);
        }
    }
    cout << "geno size: " << common_geno.size() <<","<<common_geno[0].size()<< endl;
    swap(this->myinput.geno, common_geno);
    swap(this->myinput.snpIDs, snpset);
}
void dataclient::mapping(string pheno_path, string geno_path, string pheno_pos, string geno_pos, int rowstart, int rowend, int permut)
{
    // dataclient client1;
    // initialize(sendport1, recvport1, address1, sendport2, recvport2, address2, sendport3, recvport3, address3);
    auto pre_start = chrono::high_resolution_clock::now();
    vector<vector<int64_t>> geno_scaled;
    vector<double> geno_var;
    // string pheno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/bed_template.tsv"; //gtex
    // string pheno_pos = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/phenotype/data/stableID_bed_template_v26.tsv"; // run2 gmg 
    // string pheno_pos = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/phenotype/data/stableID_bed_template_v26_reindexed.tsv"; //run2 testing
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype.tsv";
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_genotype_1kGresidualized.tsv";
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/data/genotype/re/genotype_imputed_projected_residualized_" + split_set +"_concatenated.tsv"; // scn1 and scn2 projected matrix
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run/script2/time_measure/data/projected_200_concatenated.tsv"; // For sample size testing
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/genotype/data/genotype_projected_residualized_concatenated_gmg.tsv"; //run2 scenario1 and 2 projected, residualized, concatenated
    // string geno_matrix = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/genotype/data/GOLD_concat_pca_residualized.tsv"; //run2 gold genotype for testing
    
    // string geno_pos = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/GTEx_v8_blood_WGS_variant.tsv"; //for just GTEx
    // string geno_pos = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/covariates/1KGP_data/genotype_check/remapped_1KGP_variants.tsv"; //for run2 gmg
    // string geno_pos = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/covariates/1KGP_data/genotype_check/remapped_1KGP_variants_reindexed.tsv"; //for run2 gmg
  
    this->myinput = load_genotype(pheno_pos, geno_path, geno_pos, 1000000);
    find_common(pheno_path);
    for (int r=rowstart; r<rowend; r++)
    {
        vector<string> cisVariants;
        vector<double> std_ratio;
        string geneID = this->genes[r];
        // string pheno_file = pheno_path;
        vector<double> norm_pheno = this->pheno[r];
        // read_bedfile_row(norm_pheno,geneID, pheno_file, r,4,true); // for scenario1, where is bed file
        // read_bedfile_row(norm_pheno,geneID, pheno_file, r,0,false); // for scenario2, where first column is geneID
        vector<uint64_t> range;
        vector<vector<double>> slicedgeno;
        string chromosome = this->myinput.getCisRange(geneID,range);
        int returned = this->myinput.sliceGeno(range, chromosome, -1,cisVariants, slicedgeno);
        this->owner_p1.send(returned);
        this->owner_p2.send(returned);
        this->owner_p3.send(returned);
        if (returned != 0){
            cout << string("skipping row ") + to_string(r) + string(" ") + geneID + string(" and continuing next iteration.") << endl;
            continue;
        }

        /***CENTER AND NORMALIZE VEC**/
        double row_mean = 0.0;
        double mean;
        for (int j = 0; j < norm_pheno.size(); ++j) {
            row_mean += norm_pheno[j];
        }
        this->owner_p1.send(row_mean);
        this->p1_owner.recv(mean);
        mean /= norm_pheno.size();
        double row_variance = 0.0;
        double variance;
        for (int j = 0; j < norm_pheno.size(); ++j) {
            norm_pheno[j] -= mean;
            row_variance += norm_pheno[j] * norm_pheno[j];
        }
        this->owner_p1.send(row_variance); //CP1 sums of sq_error
        this->p1_owner.recv(variance);
        double norm = sqrt(variance); 
        for (int j = 0; j < norm_pheno.size(); ++j) {
            norm_pheno[j] /= norm;
        }
        double bed_var = variance /(norm_pheno.size() - 1);
        /***CENTER AND NORMALIZE VEC**/
        vector<int64_t> pheno = ScaleVector_signed(norm_pheno, pow(10,6)); // integer version of secret CHANGE
        
        cout << string("row ") + to_string(r) + string(" ") + geneID + string(" bed file phenotype var: "+to_string(bed_var)+" and phenotype length: ") << pheno.size()<< endl;;
        /***CENTER and NORMALIZE GENO**/
        int rows = slicedgeno.size();
        int cols = slicedgeno[0].size();
        vector<vector<double>> N(rows, vector<double>(cols, 0.0));
        vector<double> row_variances(rows, 0.0);
        vector<double> geno_mean_sum(rows, 0.0); //sum of values per row
        for (int i = 0; i < rows; ++i) {
            // Calculate row mean
            double row_mean = 0.0;
            for (int j = 0; j < cols; ++j) {
                row_mean += slicedgeno[i][j];
            }
            geno_mean_sum.push_back(row_mean);
        }
        this->owner_p1.send(geno_mean_sum);
        vector<double> geno_mean;
        this->p1_owner.recv(geno_mean);
        for (int i = 0; i < rows; ++i) {
            geno_mean[i] /= cols;
            // Center and normalize the row, and calculate variance
            for (int j = 0; j < cols; ++j) {
                N[i][j] = slicedgeno[i][j] - geno_mean[i];
                row_variances[i] += N[i][j] * N[i][j];
            }
        }
        this->owner_p1.send(row_variances);
        vector<double> geno_error_sum;
        this->p1_owner.recv(geno_error_sum); //sum of squared error per row (aggregated)
        for (int i = 0; i < rows; ++i) {
            double norm = sqrt(geno_error_sum[i]);
            for (int j = 0; j < cols; ++j) {
                N[i][j] /= norm;
            }
            geno_error_sum[i] /= (cols-1); // Not the same
        }
        swap(N, slicedgeno);
        swap(geno_var, geno_error_sum);

        geno_scaled = ScaleVector(slicedgeno, pow(10,6));
 
        for (size_t j = 0; j < geno_var.size(); ++j) {
            std_ratio.push_back(sqrt(bed_var / geno_var[j]));
        }

        uint64_t inv = PowerMod(3, -1, p);
        vector<vector<uint64_t>> shares;
        for (int i = 0; i < 3; i++) {
            shares.push_back(vector<uint64_t>());
        }
        
        // making the shares: for each bit, loop through secrets
        for (int i=0; i<pheno.size(); i++)
        {
            ZZ_p x1 = random_ZZ_p();
            ZZ_p x2 = random_ZZ_p();
            ZZ_p x3;
            if (pheno[i]<0)
                x3 = (conv<ZZ_p>(pheno[i]+p)) - x1 - x2;
            else
                x3 = conv<ZZ_p>(pheno[i]) - x1 - x2;

            uint64_t send_x1 = conv<uint64_t>(x1);
            uint64_t send_x2 = conv<uint64_t>(x2);
            uint64_t send_x3 = conv<uint64_t>(x3);
            shares[0].push_back(send_x1);
            shares[0].push_back(send_x2);
            shares[1].push_back(send_x2);
            shares[1].push_back(send_x3);                    
            shares[2].push_back(send_x3);
            shares[2].push_back(send_x1);
        }
        this->owner_p1.send(shares[0]);
        this->owner_p2.send(shares[1]);
        this->owner_p3.send(shares[2]);
        
        cout << "Sent secret shared pheno values to parties.\n";
        vector<vector<uint64_t>> genoshares;
        for (int i = 0; i < 3; i++) {
            genoshares.push_back(vector<uint64_t>());
        }
        // cout << "2\n";
        // sending geno shares
        for (int i=0; i< geno_scaled.size(); i++)
        {
            for (int j=0; j< geno_scaled[0].size(); j++)
            {
                ZZ_p i1 = random_ZZ_p();
                ZZ_p i2 = random_ZZ_p();
                ZZ_p i3;
                if (geno_scaled[i][j]<0)
                    i3 = (conv<ZZ_p>(geno_scaled[i][j]+p)) - i1 - i2;
                else
                    i3 = conv<ZZ_p>(geno_scaled[i][j]) - i1 - i2;

                uint64_t send_i1 = conv<uint64_t>(i1);
                uint64_t send_i2 = conv<uint64_t>(i2);
                uint64_t send_i3 = conv<uint64_t>(i3);
                genoshares[0].push_back(send_i1);
                genoshares[0].push_back(send_i2);
                genoshares[1].push_back(send_i2);
                genoshares[1].push_back(send_i3);
                genoshares[2].push_back(send_i3);
                genoshares[2].push_back(send_i1);
            }
        }
        // cout << "1\n";
        uint64_t testg_row = geno_scaled.size();
        uint64_t testg_col = geno_scaled[0].size();
        vector<uint64_t> genoshape = {
            static_cast<uint64_t>(testg_row),
            static_cast<uint64_t>(testg_col),
        };
        this->owner_p1.send(genoshape);
        this->owner_p2.send(genoshape);
        this->owner_p3.send(genoshape);
        this->owner_p1.send(genoshares[0]);
        this->owner_p2.send(genoshares[1]);
        this->owner_p3.send(genoshares[2]);
        cout << "Sent secret shared geno values to parties.\n";
        // sending std_ratio in plaintext to party 1
        this->owner_p1.send(std_ratio);
        string serializedvariants=string(geneID+";");
        for (const std::string& str : cisVariants) {
            serializedvariants += str + ";"; // Use a suitable delimiter
        }
        int estimatedSize = static_cast<int>(serializedvariants.size());
        int bufferSize = estimatedSize;
        // Send the serialized data over the channel
        this->owner_p1.send(bufferSize);
        this->owner_p1.send(serializedvariants.data(), serializedvariants.size());
        // cout << "Geno shared " << std::time(nullptr) << endl;
        cout << "Sent secret shared geno to parties.\n";
        int complete;
        this->p1_owner.recv(complete);
        if (complete != 1)
        {
            cout << "Didn't receive confirmation. Did row finish?" << endl;
        }
    }
}

void bitdecompose(vector<uint64_t> &secrets, vector<BitVector> &bitInput)
{
    // vector<BitVector> bitInput;
    for (int i = 0; i < secrets.size(); i++)
    {
        string sinput = bitset<64>(secrets[i]).to_string();
        BitVector decomposed(sinput);
        bitInput.push_back(decomposed);
    }
}

vector<double> get_quantiles(vector<vector<double>>& phen_matrix, vector<vector<size_t>>& rank_matrix) {
    for (size_t i = 0; i < phen_matrix[0].size(); i++) {
        vector<pair<double, size_t>> sorted_indices;
        for (size_t j = 0; j < phen_matrix.size(); j++) {
            sorted_indices.emplace_back(phen_matrix[j][i], j);
        }
        sort(sorted_indices.begin(), sorted_indices.end());
        for (size_t j = 0; j < phen_matrix.size(); j++) {
            rank_matrix[i][j] = sorted_indices[j].second;
        }
    }

    size_t m = phen_matrix.size();
    size_t n = phen_matrix[0].size();

    vector<double> quantiles(m, 0.0);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            quantiles[j] += phen_matrix[rank_matrix[i][j]][i];
        }
    }

    return quantiles;
}

void sample_QN(vector<vector<double>>& phen_matrix, vector<vector<size_t>>& rank_matrix, vector<double>& total_quantiles) {
    size_t m = phen_matrix.size();
    size_t n = phen_matrix[0].size();

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            phen_matrix[rank_matrix[i][j]][i] = total_quantiles[j]/n;
        }
    }
}

vector<vector<double>> deseq2_cpm(vector<vector<uint64_t>>& counts_df) {
    int numGenes = counts_df.size();
    int numSamples = counts_df[0].size();

    vector<double> colSums(numSamples, 0.0);
    for (int j = 0; j < numSamples; ++j) {
        for (int i = 0; i < numGenes; ++i) {
            colSums[j] += static_cast<double>(counts_df[i][j]);
        }
    }
    vector<vector<double>> cpm_df(numGenes, vector<double>(numSamples, 0.0));
    for (int j = 0; j < numSamples; ++j) {
        for (int i = 0; i < numGenes; ++i) {
            cpm_df[i][j] = counts_df[i][j] / colSums[j] * 1e6;
        }
    }
    return cpm_df;
}

void dataclient::phenopreprocess(string norm_method, string pheno_input, int rowstart, int rowend, string zscorefile, vector<vector<double>>& resultVec, vector<string>& gene_string)
{
    // initialize(sendport1, recvport1, address1, sendport2, recvport2, address2, sendport3, recvport3, address3);
    auto pre_start = chrono::high_resolution_clock::now();
    vector<string> geneID;
    vector<vector<double>> pheno;
    if (norm_method == "qn")
    {
        cout << "Quantile Normalization executing... " << flush;
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_tpm_matched_filtered.tsv"; //original
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_tpm_matched_filtered_" + split_set + ".tsv";
        pheno = getTPMFromMatrixFile(pheno_input, geneID, true);
        // print_vector(testss);
        vector<vector<size_t>> rank(pheno[0].size(), vector<size_t>(pheno.size()));
        vector<double> quantiles =get_quantiles(pheno, rank);
        sample_QN(pheno, rank, quantiles);
        auto pre_end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = pre_end - pre_start;
        double durationInSeconds = duration.count();
        double durationInminutes = durationInSeconds/60.0;
        cout << durationInminutes << " minutes" << endl;
    }
    else if (norm_method == "deseq2")
    {
        cout << "Deseq2 normalization executing... ";
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_reads_matched_filtered.tsv"; // original order
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_reads_matched_filtered_" + split_set + ".tsv"; // for just GTEx, skip 2 columns
        // string matched = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/phenotype/data/all_match_filtered_read_counts_700.tsv"; // for GTEx, geuvadis, mage
        vector<vector<uint64_t>> testp = getCountFromMatrixFile(pheno_input,geneID,1);
        vector<vector<double>> cpm_df = deseq2_cpm(testp); ///CHANGE TO testp
        vector<vector<uint64_t>> phenoshares;
        for (int i = 0; i < 3; i++) {
            phenoshares.push_back(vector<uint64_t>());
        }
        uint64_t testp_row = cpm_df.size();
        uint64_t testp_col = cpm_df[0].size();
        vector<int> exclude;
        for (int i=0; i<cpm_df.size(); i++)
        {
            vector<uint64_t> share1_row, share2_row, share3_row;
            for (int j=0; j< cpm_df[0].size(); j++)
            {
                if (cpm_df[i][j]<=0)
                {
                    share1_row.clear();
                    share2_row.clear();
                    share3_row.clear();
                    testp_row--;
                    exclude.push_back(i);
                    break;
                }
                double logcount = log(cpm_df[i][j]);
                int64_t scaledcount = logcount*pow(10,5);
                // cout << scaledcount << endl;
                ZZ_p i1 = random_ZZ_p();
                ZZ_p i2 = random_ZZ_p();
                ZZ_p i3;
                if (scaledcount >= 0)
                    i3 = conv<ZZ_p>(scaledcount) - i1 - i2;
                else
                    i3 = conv<ZZ_p>(scaledcount+p) - i1 - i2;
                uint64_t send_i1 = conv<uint64_t>(i1);
                uint64_t send_i2 = conv<uint64_t>(i2);
                uint64_t send_i3 = conv<uint64_t>(i3);
                share1_row.push_back(send_i1);
                share1_row.push_back(send_i2);
                share2_row.push_back(send_i2);
                share2_row.push_back(send_i3);
                share3_row.push_back(send_i3);
                share3_row.push_back(send_i1);
            }
            phenoshares[0].insert(phenoshares[0].end(), share1_row.begin(), share1_row.end());
            phenoshares[1].insert(phenoshares[1].end(), share2_row.begin(), share2_row.end());
            phenoshares[2].insert(phenoshares[2].end(), share3_row.begin(), share3_row.end());
        }
        vector<uint64_t> matshape = {
            static_cast<uint64_t>(testp_row),
            static_cast<uint64_t>(testp_col)
        };

        this->owner_p1.send(matshape);
        this->owner_p2.send(matshape);
        this->owner_p3.send(matshape);
        this->owner_p1.send(phenoshares[0]);
        this->owner_p2.send(phenoshares[1]);
        this->owner_p3.send(phenoshares[2]);
        vector<double> ref1, ref2, ref3;
        vector<vector<double>> finalratio(testp[0].size(), vector<double>(testp.size()-exclude.size()));
        this->p1_owner.recv(ref1);
        this->p2_owner.recv(ref2);
        this->p3_owner.recv(ref3);
        int skipped = 0;
        for (int i=0; i<cpm_df.size(); i++)
        {
            auto it = find(exclude.begin(), exclude.end(), i);
            if (it != exclude.end())
            {
                // cout <<string("Gene "+to_string(i)+" excluded;\n");
                skipped++;
                continue;
            }
            for(int j=0; j<cpm_df[0].size();j++)
            {
                if (ref1[i-skipped] >=0 )
                    finalratio[j][i-skipped] = log(cpm_df[i][j]) - ref1[i-skipped]/pow(10,5);
                else
                    finalratio[j][i-skipped] = log(cpm_df[i][j]) - (ref1[i-skipped]-p)/pow(10,5);
            }
        }
        vector<double> sizefactors;
        for (int j=0; j<finalratio.size(); j++)
        {
            sort(finalratio[j].begin(), finalratio[j].end());
            size_t size = finalratio[j].size();
            if (size % 2 == 0) {
                double medval = (finalratio[j][size / 2 - 1] + finalratio[j][size / 2]) / 2.0;
                sizefactors.push_back(exp(medval));
            } else {
                sizefactors.push_back(exp(finalratio[j][size / 2]));
            }
        }
        for (int i=0; i<cpm_df.size(); i++)
        {
            for (int j=0; j<cpm_df[0].size(); j++)
            {
                cpm_df[i][j] = cpm_df[i][j]/sizefactors[j];
            }
        }
        swap(cpm_df, pheno);
        auto pre_end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = pre_end - pre_start;
        double durationInSeconds = duration.count();
        double durationInminutes = durationInSeconds/60.0;
        cout << durationInminutes << " minutes" << endl;
    }
    else 
    {
        throw invalid_argument("Please choose normalization method between qn and deseq2.\n");
    }

    for (int r=rowstart; r<rowend; r++)
    {
        if (r >= pheno.size())
        {
            cout << r << ", out of " << pheno.size() << "total lines" << endl;
            break;
        }
            
        auto rowstart = chrono::high_resolution_clock::now();
        vector<string> cisVariants;
        vector<double> std_ratio;
        int ready1, ready2, ready3;
        this->p1_owner.recv(ready1);
        this->p2_owner.recv(ready2);
        this->p3_owner.recv(ready3);
        
        if ((ready1 == 1) && (ready2 ==1) && (ready3==1))
        {
            // cout << "Client will send one phenotype\n";
            vector<BitVector> bitInput;

            
            vector<uint64_t> secrets = ScaleVector(pheno[r], pow(10,7)); // integer version of secret

            bitdecompose(secrets, bitInput);

            
            /// ZSCORE FILE
            // string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/"+zscorefile+".txt");
            // string original = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/zscores.txt");
            vector<double> zscore_input = CSVtoVector(zscorefile);
            // vector<double> pre_centered = CSVtoVector(original);
            // double pheno_var =doublevariance(pre_centered, doublemean(pre_centered));//0.9782648500530864;// 0.9596748533543238; //TODO::Don't put manual numbers here 

            auto min = min_element(zscore_input.begin(), zscore_input.end()); // CHANGE to zscore_input
            double shiftsize = abs(*min);
            // cout << "shift size: " << shiftsize << endl;
            vector<int64_t> zscores_scaled = ScaleVector_signed(zscore_input, pow(10,7));// CHANGE to zscore_input
            
            uint64_t inv = PowerMod(3, -1, p);
            vector<uint64_t> vectorsize{(uint64_t) bitInput.size(), (uint64_t) bitInput[0].size(), p, inv};
            this->owner_p1.send(vectorsize);
            this->owner_p2.send(vectorsize);
            this->owner_p3.send(vectorsize);

            vector<vector<uint64_t>> shares;
            for (int i = 0; i < 3; i++) {
                shares.push_back(vector<uint64_t>());
            }
            // making the shares: for each bit, loop through secrets
            for (int j=bitInput[0].size()-1; j >=0; j--)
            {
                for (int i=0; i<bitInput.size(); i++)
                {
                    ZZ_p x1 = random_ZZ_p();
                    ZZ_p x2 = random_ZZ_p();
                    ZZ_p x3 = to_ZZ_p(bitInput[i][j]) - x1 - x2;
                    uint64_t send_x1 = conv<uint64_t>(x1);
                    uint64_t send_x2 = conv<uint64_t>(x2);
                    uint64_t send_x3 = conv<uint64_t>(x3);
                    shares[0].push_back(send_x1);
                    shares[0].push_back(send_x2);
                    shares[1].push_back(send_x2);
                    shares[1].push_back(send_x3);                    
                    shares[2].push_back(send_x3);
                    shares[2].push_back(send_x1);
                }
            }
            this->owner_p1.send(shares[0]);
            this->owner_p2.send(shares[1]);
            this->owner_p3.send(shares[2]);
            int prelim1, prelim2, prelim3;
            this->p1_owner.recv(prelim1);
            this->p2_owner.recv(prelim2);
            this->p3_owner.recv(prelim3);
            if ((prelim1 == 1) && (prelim1 ==1) && (prelim1==1))
            {
                // cout << "Client will share zscore and identity.\n";

                vector<vector<uint64_t>> identity_shares;
                for (int i = 0; i < 3; i++) {
                    identity_shares.push_back(vector<uint64_t>());
                }
                //sending identity shares
                for (int j=0; j< bitInput.size(); j++)
                {
                    ZZ_p i1 = random_ZZ_p();
                    ZZ_p i2 = random_ZZ_p();
                    ZZ_p i3 = conv<ZZ_p>(j+1) - i1 - i2;
                    uint64_t send_i1 = conv<uint64_t>(i1);
                    uint64_t send_i2 = conv<uint64_t>(i2);
                    uint64_t send_i3 = conv<uint64_t>(i3);
                    identity_shares[0].push_back(send_i1);
                    identity_shares[0].push_back(send_i2);
                    identity_shares[1].push_back(send_i2);
                    identity_shares[1].push_back(send_i3);
                    identity_shares[2].push_back(send_i3);
                    identity_shares[2].push_back(send_i1);
                }
                this->owner_p1.send(identity_shares[0]);
                this->owner_p2.send(identity_shares[1]);
                this->owner_p3.send(identity_shares[2]);
                this->owner_p1.send(shiftsize);
                this->owner_p2.send(shiftsize);
                this->owner_p3.send(shiftsize);
                vector<vector<uint64_t>> zscore_shares;
                for (int i = 0; i < 3; i++) {
                    zscore_shares.push_back(vector<uint64_t>());
                }
                for (int i=0; i<zscores_scaled.size(); i++)
                {
                    ZZ_p z1 = random_ZZ_p();
                    ZZ_p z2 = random_ZZ_p();
                    ZZ_p z3;
                    if (zscores_scaled[i]<0)
                        z3 = (conv<ZZ_p>(zscores_scaled[i]+p)) - z1 - z2;
                    else
                        z3 = conv<ZZ_p>(zscores_scaled[i]) - z1 - z2;
                    // ZZ_p z3 = to_ZZ_p(zscores_shifted[i]) - z1 - z2;
                    uint64_t send_z1 = conv<uint64_t>(z1);
                    uint64_t send_z2 = conv<uint64_t>(z2);
                    uint64_t send_z3 = conv<uint64_t>(z3);
                    zscore_shares[0].push_back(send_z1);
                    zscore_shares[0].push_back(send_z2);
                    zscore_shares[1].push_back(send_z2);
                    zscore_shares[1].push_back(send_z3);
                    zscore_shares[2].push_back(send_z3);
                    zscore_shares[2].push_back(send_z1);
                }
                this->owner_p1.send(zscore_shares[0]);
                this->owner_p2.send(zscore_shares[1]);
                this->owner_p3.send(zscore_shares[2]);
            }
            
            // cout << "Secrets shared " << std::time(nullptr) << endl;
            // cout << "Sent secret shared row values to parties.\n";
        }
       
        vector<uint64_t> invcdf_1, invcdf_2, invcdf_3;
        vector<double> row_result;
        this->p1_owner.recv(invcdf_1);
        this->p2_owner.recv(invcdf_2);
        this->p3_owner.recv(invcdf_3);
        vector<ZZ_p> result1 = convVec(invcdf_1);
        vector<ZZ_p> result2 = convVec(invcdf_2);
        vector<ZZ_p> result3 = convVec(invcdf_3);
        for (int i=0;i<invcdf_1.size()/2; i++)
        {
            int64_t unshifted;
            ZZ_p final = result1[2*i]+result2[2*i]+result3[2*i];
            if (conv<uint64_t>(final) > p/2)
                unshifted = conv<int64_t>(final) - p;
            else
                unshifted = conv<int64_t>(final);
            double final_res = static_cast<double>(unshifted)/pow(10,7);
            row_result.push_back(final_res);
        }
        resultVec.push_back(row_result);
        gene_string.push_back(geneID[r]);
        auto rowend = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = rowend - rowstart;
        double durationInSeconds = duration.count();
        double durationInminutes = durationInSeconds/60.0;
        cout << "Row "<< r << " execution time: " << durationInminutes << " minutes" << endl;
    }
}
void dataclient::close()
{
    std::cout <<"Dataclient closing" << std::endl;
    this->owner_p1.close();
    this->owner_p2.close();
    this->owner_p3.close();
    this->p1_owner.close();
    this->p2_owner.close();
    this->p3_owner.close();
    std::cout << "Channels closed." << std::endl;
}
void dataclient_mapping(string pheno_path, string geno_path, string pheno_pos, string geno_pos, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend, int permut)
{
    dataclient client1;
    client1.initialize(sendport1, recvport1, address1, sendport2, recvport2, address2, sendport3, recvport3, address3);
    // client1.find_common(pheno_path);
    client1.mapping(pheno_path, geno_path, pheno_pos, geno_pos, rowstart, rowend, permut);
    client1.close();
}
void dataclient_preprocessing(string norm_method, string pheno_input, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,string zscorefile, vector<vector<double>>& resultVec, vector<string>& gene_string)
{
    cout << "owner calling privateQTLII script" << endl;
    dataclient client2;
    client2.initialize(sendport1, recvport1, address1, sendport2, recvport2, address2, sendport3, recvport3, address3);
    client2.phenopreprocess(norm_method, pheno_input, rowstart, rowend, zscorefile, resultVec, gene_string);
    client2.close();
}