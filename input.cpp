#include "input.h"
#include <typeinfo>
#include <chrono>

// Data owners convert plink bed file to dataframe using gtex pipeline
void getGenotype(const string& filename, vector<vector<double>>& geno, vector<string>& snpID) {
    trace("getGenotype");
    // cout << "Fetching genotype into memory..\n";
    vector<vector<double>> rowsData;
    ifstream data(filename);
    string line;
    // int currentRow = 0;

    // Skip the first row (header)
    getline(data, line);

    while (getline(data, line)) {
        
        stringstream lineStream(line);
        string cell;
        getline(lineStream, cell, '\t');
        snpID.push_back(cell);
        // getline(lineStream, cell, '\t');

        vector<double> rowVector;
        while (getline(lineStream, cell, '\t')) {
            try {
                double entry = stod(cell);
                rowVector.push_back(entry);
            } catch (const exception& e) {
                cout << "Error in getgenotype" << endl;
                cerr << "Exception caught: " << e.what() << endl;
            }
        }
        rowsData.push_back(rowVector);
    }

    // Close the file after reading
    data.close();

    swap(geno, rowsData);
}
 
// gene positions
unordered_map<string, GeneData> readDataFile(string& filename) {
    trace("readDataFile");
    unordered_map<string, GeneData> geneMap;
    ifstream data(filename);
    string line;
    // Skip the header row
    getline(data, line);

    while (getline(data, line)) {
        istringstream lineStream(line);
        string geneID, chr;
        uint64_t start, end;
        lineStream >> geneID>> chr;
        if (chr.size() >= 3 && chr.substr(0, 3) == "chr") {
            chr = chr.substr(3); // Strip "chr" from the beginning
        }
        lineStream >> start >> end;
        GeneData geneData;
        geneData.chr = chr;
        geneData.start = start;
        geneData.end = end;
        geneMap[geneID] = geneData;
    }

    data.close();
    return geneMap;
}

void getSNPpos(string& geno_pos, vector<uint64_t>& positions, vector<string>& chrpos)
{
    trace("getSNPpos");
    ifstream data(geno_pos);
    string line;
    // Skip the header row
    getline(data, line);
    while (getline(data, line)) 
    {
        istringstream lineStream(line);
        string snpID, chr;
        uint64_t pos;
        lineStream >> snpID >> chr >> pos;
        positions.push_back(pos); //ordered positions of SNPs
        chrpos.push_back(chr);
    }
}
prepareInput::prepareInput(string& pheno_pos, string& geno_matrix, string& geno_pos, uint64_t window)
{
    trace("prepareInput");
    ciswindow = window;
    cout << "loading data..." << endl;
    genePos = readDataFile(pheno_pos);
    cout << "gene positions loaded" << endl;
    getSNPpos(geno_pos, snpPos, snpChr);
    cout << "snp positions loaded" << endl;
    getGenotype(geno_matrix, geno, snpIDs);
    cout << string("genotype matrix loaded: "+to_string(geno.size())+","+to_string(geno[0].size())+"\n");
}
string prepareInput::getCisRange(string geneID, vector<uint64_t>& positions)
{ 
    trace("getCisRange");
    cout << "getting cis range:";
    // get cis range based on gene position
    const GeneData& geneData = genePos[geneID];
    positions.push_back(geneData.start);
    positions.push_back(geneData.end);
    // positions = {geneData.start, geneData.end};
    cout << " done." << endl;
    return geneData.chr;
}
vector<uint64_t> prepareInput::getSNPrange(uint64_t start, uint64_t end, string chrnum, vector<string>& cisSnpIds)
{     
    trace("prepareInput::getSNPrange");
    // get snp indices that fall within range
    vector<uint64_t> indices;
    auto lb_it = lower_bound(snpPos.begin(), snpPos.end(), static_cast<int64_t>(start) - static_cast<int64_t>(ciswindow));
    auto ub_it = upper_bound(snpPos.begin(), snpPos.end(), end + ciswindow);

    // Get the actual lower bound value in snpPos vector
    uint64_t lower_bound_value = (lb_it != snpPos.end()) ? *lb_it : *min_element(snpPos.begin(), snpPos.end());

    // Get the actual upper bound value in snpPos vector
    uint64_t upper_bound_value = (ub_it != snpPos.end()) ? *ub_it : *max_element(snpPos.begin(), snpPos.end());
    for (size_t i = 0; i < geno.size(); ++i) {
        int64_t cisdistance = snpPos[i] - start;
        if(abs(cisdistance)<=ciswindow && snpChr[i] == chrnum)
        {
            indices.push_back(i);
            cisSnpIds.push_back(snpIDs[i]);
        }
    }

    cout << "done." << endl;
    return indices;
}
int prepareInput::sliceGeno(vector<uint64_t> positions, string& chr, int64_t missing, vector<string>& cisSNPs, vector<vector<double>>& slicedmatrix)
{
    trace("prepareInput::sliceGeno");
    cout << "in slicegeno...";
    uint64_t start = positions[0];
    uint64_t end = positions[1];
    vector<uint64_t> idx = getSNPrange(start, end, chr, cisSNPs);
    // vector<vector<double>> slicedmatrix;
    // cout << "idx size " << idx.size() << endl;
    if (idx.size() == 0)
    {
        cout << "\nthis gene has no SNPs in cis range." << endl;
        return 1;
    }
    else
    {
        for (uint64_t index : idx)
        {   
            int sum = 0;
            vector<int> missing_idx;
            for (int j=0; j < geno[index].size(); j++)
            {
                if (geno[index][j]==missing)
                {
                    missing_idx.push_back(j);
                    continue;
                }
                sum+=geno[index][j];
            }
            // cout << "missing_index " << missing_idx.size() << endl;
            double avg = static_cast<double>(sum)/(geno[index].size()-missing_idx.size());
            for (int i=0; i<missing_idx.size(); i++)
            {
                geno[index][missing_idx[i]] = avg;
            }
            slicedmatrix.push_back(geno[index]);
        }
        cout << "slicing complete\n";
        return 0;
    }
    // int head = 0;
    // cout << "geno size: " << geno.size() << endl;
    
}

