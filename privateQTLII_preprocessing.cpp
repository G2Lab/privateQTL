
#include "privateQTL_cp.h"
#include <cmath>
#include "utils.h"
#include <chrono>
#include "input.h"
#include <sys/resource.h>
struct rusage r_usage;

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

void dataclient(string norm_method, string split_set, int sendport1, int recvport1, string address1, int sendport2, int recvport2, string address2, int sendport3, int recvport3, string address3,int rowstart, int rowend,string zscorefile, vector<vector<double>>& resultVec, vector<string>& gene_string)
{
    IOService ios;
    Endpoint epsend1(ios, address1, sendport1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendport2, EpMode::Server);
    Endpoint epsend3(ios, address3, sendport3, EpMode::Server);
    Endpoint eprecv1(ios, address1, recvport1, EpMode::Client);
    Endpoint eprecv2(ios, address2, recvport2, EpMode::Client);
    Endpoint eprecv3(ios, address3, recvport3, EpMode::Client);
    Channel owner_p1 = epsend1.addChannel();
    Channel owner_p2 = epsend2.addChannel();
    Channel owner_p3 = epsend3.addChannel();
    Channel p1_owner = eprecv1.addChannel();
    Channel p2_owner = eprecv2.addChannel();
    Channel p3_owner = eprecv3.addChannel();
    cout << "Owner established channels with the computing parties.\n";
    NTL::SetSeed((NTL::conv<NTL::ZZ>((long)27))); // Seed change
    vector<vector<double>> pheno;
    uint64_t p = pow(2,50);
    ZZ_p::init(to_ZZ(p)); 
    cout << "P: " << p << endl;
    owner_p1.send(p);
    owner_p2.send(p);
    owner_p3.send(p);
    auto pre_start = chrono::high_resolution_clock::now();
    vector<string> geneID;
    if (norm_method == "qn")
    {
        cout << "Quantile Normalization executing... " << flush;
        // string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_tpm_matched_filtered.tsv"; //original
        string matched = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/blood_tpm_matched_filtered_" + split_set + ".tsv";
        pheno = getTPMFromMatrixFile(matched, geneID);
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
        string matched = "/gpfs/commons/groups/gursoy_lab/ykim/QTL_proj/run2/phenotype/data/all_match_filtered_read_counts_700.tsv"; // for GTEx, geuvadis, mage
        vector<vector<uint64_t>> testp = getCountFromMatrixFile(matched,geneID,1);

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

        owner_p1.send(matshape);
        owner_p2.send(matshape);
        owner_p3.send(matshape);

        owner_p1.send(phenoshares[0]);
        owner_p2.send(phenoshares[1]);
        owner_p3.send(phenoshares[2]);
        vector<double> ref1, ref2, ref3;
        vector<vector<double>> finalratio(testp[0].size(), vector<double>(testp.size()-exclude.size()));
        p1_owner.recv(ref1);
        p2_owner.recv(ref2);
        p3_owner.recv(ref3);
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

    vector<vector<int64_t>> geno_scaled;
    vector<double> geno_var, pheno_var;
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
        p1_owner.recv(ready1);
        p2_owner.recv(ready2);
        p3_owner.recv(ready3);
        
        if ((ready1 == 1) && (ready2 ==1) && (ready3==1))
        {
            cout << "Client will send one phenotype\n";
            vector<BitVector> bitInput;

            
            vector<uint64_t> secrets = ScaleVector(pheno[r], pow(10,7)); // integer version of secret

            bitdecompose(secrets, bitInput);

            
            /// ZSCORE FILE
            string zscore_filename = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/"+zscorefile+".txt");
            string original = string("/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/rnaseq/data/blood/zscores.txt");
            vector<double> zscore_input = CSVtoVector(zscore_filename);
            vector<double> pre_centered = CSVtoVector(original);
            double pheno_var =doublevariance(pre_centered, doublemean(pre_centered));//0.9782648500530864;// 0.9596748533543238; //TODO::Don't put manual numbers here 

            auto min = min_element(zscore_input.begin(), zscore_input.end()); // CHANGE to zscore_input
            double shiftsize = abs(*min);
            // cout << "shift size: " << shiftsize << endl;
            vector<int64_t> zscores_scaled = ScaleVector_signed(zscore_input, pow(10,7));// CHANGE to zscore_input
            
            uint64_t inv = PowerMod(3, -1, p);
            vector<uint64_t> vectorsize{(uint64_t) bitInput.size(), (uint64_t) bitInput[0].size(), p, inv};
            owner_p1.send(vectorsize);
            owner_p2.send(vectorsize);
            owner_p3.send(vectorsize);

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
            owner_p1.send(shares[0]);
            owner_p2.send(shares[1]);
            owner_p3.send(shares[2]);
            int prelim1, prelim2, prelim3;
            p1_owner.recv(prelim1);
            p2_owner.recv(prelim2);
            p3_owner.recv(prelim3);
            if ((prelim1 == 1) && (prelim1 ==1) && (prelim1==1))
            {
                cout << "Client will share zscore and identity.\n";

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
                owner_p1.send(identity_shares[0]);
                owner_p2.send(identity_shares[1]);
                owner_p3.send(identity_shares[2]);
                owner_p1.send(shiftsize);
                owner_p2.send(shiftsize);
                owner_p3.send(shiftsize);
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
                owner_p1.send(zscore_shares[0]);
                owner_p2.send(zscore_shares[1]);
                owner_p3.send(zscore_shares[2]);
            }
            
            // cout << "Secrets shared " << std::time(nullptr) << endl;
            // cout << "Sent secret shared row values to parties.\n";
        }
       
        vector<uint64_t> invcdf_1, invcdf_2, invcdf_3;
        vector<double> row_result;
        p1_owner.recv(invcdf_1);
        p2_owner.recv(invcdf_2);
        p3_owner.recv(invcdf_3);
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

    owner_p1.close();
    owner_p2.close();
    owner_p3.close();
    p1_owner.close();
    p2_owner.close();
    p3_owner.close();
    epsend1.stop();
    epsend2.stop();
    epsend3.stop();
    eprecv1.stop();
    eprecv2.stop();
    eprecv3.stop();
    ios.stop();
}

vector<ZZ_p> convert(vector<int> input)
{
    vector<ZZ_p> output;
    for (int i=0; i<input.size(); i++)
    {
        output.push_back(conv<ZZ_p>(input[i]));
    }
    return output;
}
void runMPC(string norm_method, int pid,  string ownerIP, int ownerPort, int toOwnerPort,  string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2,int rowstart, int rowend)
{    
    mpc CP_i;
    vector<vector<double>> resultVectors;
    CP_i.initialize(pid, ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2);
    if (norm_method == "deseq2")
    {
        CP_i.receivePheno();
        CP_i.logRatio();
    }
    auto invcdf_start = chrono::high_resolution_clock::now();
    for (int row=rowstart; row<rowend; row++)
    {
        CP_i.ready();
        CP_i.receiveSecrets();
        CP_i.genperm(row, norm_method);
        CP_i.clearVectors();
    }
    auto invcdf_end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = invcdf_end - invcdf_start;
    double durationInSeconds = duration.count();
    double durationInminutes = durationInSeconds/60.0;
    CP_i.close();
}
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
void startMPCset(string norm_method, string split_set, vector<int>& availPorts, int startidx, vector<string>& address, int startrow, int endrow, string zscorefile, int CPU_core, vector<vector<double>>& resultVec, vector<string>& gene_string)
{
    vector<thread> threads;
    thread dataClientThread([=,  &resultVec, &gene_string]() {
        dataclient(norm_method, split_set, availPorts[startidx + 6], availPorts[startidx + 9], address[1], availPorts[startidx + 7], availPorts[startidx + 10], address[2], 
        availPorts[startidx + 8], availPorts[startidx + 11], address[3], startrow, endrow, zscorefile, resultVec, gene_string);
    });
    setThreadAffinity(dataClientThread,CPU_core+0);
    threads.emplace_back(move(dataClientThread));

    thread runMPC1([=]() {
        runMPC(norm_method, 0, address[0], availPorts[startidx + 6], availPorts[startidx + 9], address[2], availPorts[startidx + 0], availPorts[startidx + 2], address[3], availPorts[startidx + 1], availPorts[startidx + 4], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC1,CPU_core+1);
    threads.emplace_back(move(runMPC1));
    thread runMPC2([=]() {
        runMPC(norm_method, 1, address[0], availPorts[startidx + 7], availPorts[startidx + 10], address[3], availPorts[startidx + 3], availPorts[startidx + 5], address[1], availPorts[startidx + 2], availPorts[startidx + 0], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC2,CPU_core+2);
    threads.emplace_back(move(runMPC2));

    thread runMPC3([=]() {
        runMPC(norm_method, 2, address[0], availPorts[startidx + 8], availPorts[startidx + 11], address[1], availPorts[startidx + 4], availPorts[startidx + 1], address[2], availPorts[startidx + 5], availPorts[startidx + 3], 
        startrow, endrow);
    });
    setThreadAffinity(runMPC3,CPU_core+2);
    threads.emplace_back(move(runMPC3));
    for (auto& thread : threads) {
        thread.join();
    }
}

int main(int argc, char* argv[])
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
    string zscorefile = argv[4];
    string norm_method = argv[5];
    string split_set = argv[6];


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
        startMPCset(norm_method, split_set, openPorts, 0, address, lo_row, mid_row, zscorefile, 0,resultVec1,gene_string1);
    });
    thread thread2([&]() {
        startMPCset(norm_method, split_set, openPorts, 12, address, mid_row, hi_row,zscorefile, 0,resultVec2,gene_string2);
    });

    threads.emplace_back(move(thread1));
    threads.emplace_back(move(thread2));

    for (auto& thread : threads) {
        thread.join();
    }
    resultVec1.insert(resultVec1.end(), resultVec2.begin(), resultVec2.end());
    gene_string1.insert(gene_string1.end(), gene_string2.begin(), gene_string2.end());
    writeNormalizedToTSV(resultVec1, gene_string1, "private_deseq2_invcdf_gmg_normalized");
    
    int siteA_n;
    int siteB_n;
    int siteC_n;
    if (split_set == "set1") ///CHANGE THIS
    {
        siteA_n = 250;//300;
        siteB_n = 300;//250;
        siteC_n = 150;//120;
    }
    else
    {
        siteA_n = 300;
        siteB_n = 300;
        siteC_n = 70;
    }
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
    writematrixToTSV(siteA_cov, "private_deseq2_invcdf_gmg_siteA_pc");
    writematrixToTSV(siteB_cov, "private_deseq2_invcdf_gmg_siteB_pc");
    writematrixToTSV(siteC_cov, "private_deseq2_invcdf_gmg_siteC_pc");
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
    string to_write = "private_deseq2_invcdf_"+split_set+"_residualized_700"; //CHANGE
    writeNormalizedToTSV(res_A, gene_string1, to_write);
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