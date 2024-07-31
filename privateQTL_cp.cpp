#include "privateQTL_cp.h"

const static Eigen::IOFormat TSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

void mpc::writeEigenToTSV(Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& mat, const string& name)
{
    string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + ".tsv";
    ofstream file(filename);
    if (file.is_open()) {
        file << mat.format(TSVFormat);
        file.close();
        std::cout << name + " matrix successfully written to TSV file.\n";
    } else {
        std::cout << "Error opening the file." << std::endl;
    }
}
double cal_accuracy(vector<double>& predicted, string& filename, int N)
{
    vector<double> actual = getRowFromMatrixFile(filename,N);
    // vector<double> actual = CSVtoVector(filename);
    // cout << string("row: "+to_string(N)+", numCol: "+to_string(numCol)+"\n");
    if (predicted.size() != actual.size()) {
        cout << "predicted size: " + to_string(predicted.size()) << endl;
        cout << "actual size: " + to_string(actual.size()) << endl;
        // throw runtime_error("Vector sizes do not match");
    }

    double sumSquaredDiff = 0.0;
    for (size_t i = 0; i < predicted.size(); ++i) {
        double diff = predicted[i] - actual[i];
        sumSquaredDiff += diff * diff;
    }

    double mse = sumSquaredDiff / predicted.size();
    double rmse = sqrt(mse);
    // cout << string("RMSE: " + to_string(rmse)) << endl;
    auto max = max_element(actual.begin(), actual.end());
    auto min = min_element(actual.begin(), actual.end());
    double normalized = *max - *min;
    return 1.0 - rmse/normalized;

}
bool mpc::initialize(int pid, string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    this->pid = pid;
    globalprng.SetSeed(commonSeed);
    if (!setupChannels(ownerIP, ownerPort, toOwnerPort, address1, recPort1, sendPort1, address2, recPort2, sendPort2))
    {
        cout << "Network channel failed.\n";
    }
    if (!setupSeeds())
    {
        cout << "PRNG setup failed.\n";
    }
    this->dataowner.recv(this->p);
    ZZ_p::init(to_ZZ(this->p));
    this->p = p;
    return true;
}
bool mpc::setupChannels(string ownerIP, int ownerPort, int toOwnerPort, string address1, int recPort1, int sendPort1, string address2, int recPort2, int sendPort2)
{
    Endpoint p_owner(ios, ownerIP, ownerPort, EpMode::Client);
    Endpoint ownersend(ios, ownerIP, toOwnerPort, EpMode::Server);
    Endpoint eprec1(ios, address1, recPort1, EpMode::Client);
    Endpoint eprec2(ios, address2, recPort2, EpMode::Client);
    Endpoint epsend1(ios, address1, sendPort1, EpMode::Server);
    Endpoint epsend2(ios, address2, sendPort2, EpMode::Server);

    this->dataowner = p_owner.addChannel();
    this->toOwner = ownersend.addChannel();
    this->fromPlus = eprec1.addChannel();
    this->fromMinus = eprec2.addChannel();
    this->toPlus = epsend1.addChannel();
    this->toMinus = epsend2.addChannel();
    cout << "Established channels with the computing parties and data owner.\n";
    p_owner.stop();
    ownersend.stop();
    eprec1.stop();
    eprec2.stop();
    epsend1.stop();
    epsend2.stop();
    // ios.stop();
    return true;
}
bool mpc::setupSeeds()
{
    // cout << "setting up shared seeds within parties..\n";
    // lower number pid party sends prng seed to higher number pid party
    //  Check if they already share a prng
    auto it = this->seedpair.find((this->pid + 1) % 3);
    auto it2 = this->seedpair.find((this->pid + 2) % 3);
    if (it != this->seedpair.end())
    {
        cout << "Party " << this->pid << " already share a prng with party " << (this->pid + 1) % 3 << endl;
        return false;
    }
    if (it2 != this->seedpair.end())
    {
        cout << "Party " << this->pid << " already share a prng with party " << (this->pid + 2) % 3 << endl;
        return false;
    }
    try
    {
        this->localprng = new PRNG(oc::toBlock(this->pid)); //NEEDS CHANGE

        uint64_t plusSeed = this->localprng->get<uint64_t>();
        this->toPlus.send(plusSeed);

        uint64_t minusSeed;
        this->fromMinus.recv(minusSeed);

        PRNG *wPlus = new PRNG(oc::toBlock(plusSeed));
        PRNG *wMinus = new PRNG(oc::toBlock(minusSeed));

        /// add prng to each party map
        this->seedpair.insert({(this->pid + 1) % 3, wPlus});
        this->seedpair.insert({(this->pid + 2) % 3, wMinus});
    }
    catch (const std::exception &e)
    {
        cout << "except\n";
        cout << e.what() << std::endl;
    }

    return true;
}

void mpc::ready()
{
    bool empty = areVectorsEmpty();
    // if (empty)
        // cout << string(to_string(this->pid) + " ready to receive.\n");
    int ready = static_cast<int>(empty);
    this->toOwner.send(ready);
}
void mpc::receiveSecrets()
{
    vector<uint64_t> dest;
    this->dataowner.recv(dest);
    this->n = dest[0]; //number of secrets
    this->lk = dest[1]; //number of bits
    this->inv = dest[3]; //inv 
    vector<uint64_t> result;
    try 
    {
        // receiving shares;
        this->dataowner.recv(result);
        vector<ZZ_p> temp(result.size());
        temp = convVec(result);
        // vector<ZZ_p> temp(result.begin(), result.end());
        swap(this->shares, temp);
        // cout <<"secrets received.\n";
        bool prelim = this->identity.empty() && this->zscores.empty();
        int ready = static_cast<int>(prelim);
        this->toOwner.send(ready);
        if(prelim)
        {   
            // cout << string(to_string(this->pid) + " ready to receive identity and zscores.\n");
            //receiving identity shares
            vector<uint64_t> identity_result;
            this->dataowner.recv(identity_result);
            vector<ZZ_p> identity_temp(identity_result.begin(), identity_result.end());
            swap(this->identity, identity_temp);
            
            // receiving zscore shares
            double shiftsize;
            this->dataowner.recv(shiftsize);
            this->shiftsize=shiftsize;

            vector<uint64_t> zscore_result;
            this->dataowner.recv(zscore_result);
            vector<ZZ_p> zscore_temp(zscore_result.begin(), zscore_result.end());
            swap(this->zscores, zscore_temp);
        }
        
    }
    catch (const exception &e)
    {
        cout << "except\n";
        cout << e.what() << endl;
    }
    // cout <<this->pid <<endl;
    // print_vector(this->identity);
    // cout << string(to_string(this->pid) + " Secrets received.\n");
}

// template class std::vector<ZZ_p>;

void mpc::assertSize(vector<ZZ_p> pi, string tag) {
    for (ZZ_p z : pi) {
        uint64_t value = conv<uint64_t>(z);
        if (value > pi.size() || value < 1) {
            string output = string("[" + to_string(this->pid) + "][" + tag + "] Error! value: " + to_string(value) + " size: " + to_string(pi.size()) + "\n");
            cout << output;
            throw logic_error("["+tag+"]: permutation is bigger than vector size.");
        }
    }    
}
vector<uint64_t> mpc::apply_plaintext_perm(vector<uint64_t> rho, vector<uint64_t> sigma)
{
    //Applying rho on sigma rho o sigma
    if (rho.size() != sigma.size())
        throw invalid_argument("This is regular apply permutation; both pis should be same size.");
    vector<uint64_t> appliedvec(rho.size());
    for (int i=0; i < sigma.size(); i++)
    {
        appliedvec[rho[i]-1] = sigma[i];
    }
    return appliedvec;
}

vector<ZZ_p> mpc::Frand(uint64_t bufferSize) 
{
    if (bufferSize <= 1)
    {
        throw std::invalid_argument("bufferSize must be positive");
    }
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &plus_PRNG = *(plus_it->second);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    vector<uint64_t> pi(2*bufferSize);
    for (uint64_t i = 0; i < bufferSize; i++)
    {
        pi[i*2+0] = i + 1;
        pi[i*2+1] = i + 1;
    }
    for (uint64_t i = bufferSize - 1; i >= 1; i--)
    {
        // random index chosen to be swapped
        uint64_t j = minus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp = pi[i*2+0];
        pi[i*2+0] = pi[j*2+0];
        pi[j*2+0] = temp;

        uint64_t k = plus_PRNG.get<uint64_t>() % (i + 1);
        uint64_t temp2 = pi[i*2+1];
        pi[i*2+1] = pi[k*2+1];
        pi[k*2+1] = temp2;
    } 
    vector<ZZ_p> ZZ_pi(pi.begin(), pi.end());
    assertSize(ZZ_pi, "Frand");
    return ZZ_pi;
}

vector<ZZ_p> mpc::reveal(vector<ZZ_p>& pi, bool isperm)
{
    // cout <<"reveal entered.\n";
    vector<uint64_t> share1;
    vector<uint64_t> share2;
    for (int i=0; i<pi.size(); i+=2)
    {
        uint64_t r = conv<uint64_t>(pi[i]);
        uint64_t q = conv<uint64_t>(pi[i+1]);
        if(r >= this->p) {
            cout << "r too big: " << r << endl;
        }
        if(q >= this->p) {
            cout << "q too big: " << q << endl;
        }
        share1.push_back(r);
        share2.push_back(q);
    }
    this->toPlus.send(share1);
    this->toMinus.send(share2);
    vector<uint64_t> receivedMinus, receivedPlus;
    this->fromMinus.recv(receivedMinus);
    this->fromPlus.recv(receivedPlus);
    
    
    if (receivedMinus != receivedPlus) 
    {
        cout << "error! in " << this->pid << endl;
        print_vector(receivedMinus);
        print_vector(receivedPlus);
        try 
        {
            throw logic_error("received inputs do not match."); // <-- THIS IS SOMETIMES THROWN BUT NOT CAUGHT
        } 
        catch (const std::string& ex) {
            std::cerr << "Exception occurred: " << ex << std::endl;
        }
        close();
    }

    vector<ZZ_p> reconstructed;
    vector<uint64_t> final(share1.size());
    vector<uint64_t> temp(share1.size());
    if (isperm)
    {
        if (this->pid==0)
        {
            temp = apply_plaintext_perm(share2, share1);
            final = apply_plaintext_perm(receivedMinus, temp);
        }
        else if (this->pid==1)
        {
            temp=apply_plaintext_perm(share1, receivedMinus);
            final=apply_plaintext_perm(share2, temp);
        }
        else
        {
            temp=apply_plaintext_perm(receivedMinus,share2);
            final=apply_plaintext_perm(share1, temp);
        }
        reconstructed = convVec(final);
    }

    else 
    {
        for (int i=0; i<receivedMinus.size(); i++)
        {
            if (receivedMinus[i] >= this->p) {
                cout << "too big! Minus: " << receivedMinus[i] << endl;
            }
            reconstructed.push_back(conv<ZZ_p>(receivedMinus[i])+conv<ZZ_p>(share1[i])+conv<ZZ_p>(share2[i]));
            // reconstructed.push_back(conv<ZZ_p>(receivedMinus[i]));
            uint64_t result = conv<uint64_t>(reconstructed[i]);
            if (result >= this->p) {
                cout << "too big! Result: " << result << endl;
            }
        }
    }
    return reconstructed;
    
}

void mpc::reshare(vector<ZZ_p>& shares, int reshareID)
{
    vector<ZZ_p> randoms(3);
    vector<uint64_t> sendvec;
    if (reshareID == (this->pid + 1) % 3) // share to plus
    {
        auto minus_it = this->seedpair.find((this->pid + 2) % 3);
        PRNG &minus_PRNG = *(minus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(minus_PRNG.get<uint64_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int j=0; j<shares.size()/2; j++)
        {
            shares[2*j]+=randoms[this->pid];
            shares[2*j+1]+=randoms[(this->pid +1)%3];
            sendvec.push_back(conv<uint64_t>(shares[2*j+1]));
            // this->toPlus.send(conv<uint64_t>(shares[2*j+1]));
        }
        this->toPlus.send(sendvec);
    }
    else if (reshareID == (this->pid + 2) % 3) //share to minus
    {
        auto plus_it = this->seedpair.find((this->pid + 1) % 3);
        PRNG &plus_PRNG = *(plus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(plus_PRNG.get<uint64_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int j=0; j<shares.size()/2; j++)
        {
            shares[2*j]+=randoms[this->pid];
            shares[2*j+1]+=randoms[(this->pid +1)%3];
            sendvec.push_back(conv<uint64_t>(shares[2*j]));
        }
        this->toMinus.send(sendvec);
    }
    else 
    {
       vector<ZZ_p> newshare(shares.size());
       vector<uint64_t> minus(shares.size()/2), plus(shares.size()/2);
       this->fromMinus.recv(minus);
       this->fromPlus.recv(plus);
        for (int i=0; i<shares.size()/2; i++)
        {
            newshare[2*i] = to_ZZ_p(minus[i]);
            newshare[2*i+1] = to_ZZ_p(plus[i]);
            // this->fromMinus.recv(newshare[2*i]);
            // this->fromPlus.recv(newshare[2*i+1]);
        }
        // vector<ZZ_p> temp(newshare.begin(), newshare.end());
        shares.swap(newshare); 

    }
}
void mpc::reshareM(vector<vector<ZZ_p>>& shares, int reshareID)
{
    vector<ZZ_p> randoms(3);
    vector<vector<uint64_t>> sendvec(shares.size(), vector<uint64_t>(shares[0].size()));
    if (reshareID == (this->pid + 1) % 3) // share to plus
    {
        auto minus_it = this->seedpair.find((this->pid + 2) % 3);
        PRNG &minus_PRNG = *(minus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(minus_PRNG.get<uint64_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int i=0; i<shares.size(); i++)
        {
            for (int j=0; j<shares[0].size()/2; j++)
            {
                shares[i][2*j]+=randoms[this->pid];
                shares[i][2*j+1]+=randoms[(this->pid +1)%3];
                sendvec[i][j] = conv<uint64_t>(shares[i][2*j+1]);
            }
            this->toPlus.send(sendvec[i]); 
        }
    }
    else if (reshareID == (this->pid + 2) % 3) //share to minus
    {
        auto plus_it = this->seedpair.find((this->pid + 1) % 3);
        PRNG &plus_PRNG = *(plus_it->second);
        for (int i=0; i<2; i++)
        {
            randoms[i] = to_ZZ_p(plus_PRNG.get<uint64_t>());
        }
        randoms[2] = -randoms[0]-randoms[1];
        for (int i=0; i<shares.size(); i++)
        {
            for (int j=0; j<shares[0].size()/2; j++)
            {
                shares[i][2*j]+=randoms[this->pid];
                shares[i][2*j+1]+=randoms[(this->pid +1)%3];
                sendvec[i][j] = conv<uint64_t>(shares[i][2*j]);
            }
            this->toMinus.send(sendvec[i]);
        }
    }
    else 
    {
       vector<vector<ZZ_p>> newshare(shares.size(), vector<ZZ_p>(shares[0].size()));
        for(int i=0; i<shares.size(); i++)
        {
            vector<uint64_t> minus(shares[0].size());
            vector<uint64_t> plus(shares[0].size());
            this->fromMinus.recv(minus);
            this->fromPlus.recv(plus);
            for(int j=0; j<shares[0].size()/2; j++)
            {
                newshare[i][2*j] = to_ZZ_p(minus[j]);
                newshare[i][2*j+1] = to_ZZ_p(plus[j]);
            }
        }
        // vector<vector<ZZ_p>> temp(newshare.begin(), newshare.end());
        shares.swap(newshare); 
    }
}
//is adding alpha beta gamma necessary? 
vector<ZZ_p> mpc::Fmult(vector<ZZ_p> k_i, vector<ZZ_p> s_i)
{
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    PRNG &plus_PRNG = *(plus_it->second);
    vector<ZZ_p> t_i(2);
    ZZ_p ri = k_i[0]*s_i[0] + k_i[1]*s_i[0] + k_i[0]*s_i[1];
    ri+=to_ZZ_p(minus_PRNG.get<uint64_t>()); // result + r1 - r2;
    ri-=to_ZZ_p(plus_PRNG.get<uint64_t>());
    t_i[0]=ri;
    this->toMinus.send(conv<uint64_t>(ri));
    uint64_t received_data;
    this->fromPlus.recv(received_data);
    conv(t_i[1], received_data);
    return t_i;
}

vector<ZZ_p> mpc::genbitperm(vector<ZZ_p> &keybit)
{
    vector<ZZ_p> f0(keybit.size()), f1(keybit.size()), s0(keybit.size()), s1(keybit.size());
    vector<ZZ_p> s(2, to_ZZ_p(0));
    vector<ZZ_p> p;
    f1 = keybit;
    for (int i=0; i<keybit.size()/2; i++)
    {
        f0[2*i] = this->inv-keybit[2*i];
        f0[2*i+1] = this->inv-keybit[2*i+1];
        s[0]+=f0[2*i];
        s[1]+=f0[2*i+1];
        s0[2*i] = s[0];
        s0[2*i+1] = s[1];
    }
    for (int i=0; i<keybit.size()/2; i++)
    {
        s[0]+=f1[2*i];
        s[1]+=f1[2*i+1];
        s1[2*i] = s[0];
        s1[2*i+1] = s[1];
    }
    vector<ZZ_p> t;
    for (int i=0; i <keybit.size()/2; i++)
    {
        vector<ZZ_p> sminus(2);
        vector<ZZ_p> ti(2);
        sminus[0] = s1[2*i] - s0[2*i];
        sminus[1] = s1[2*i+1] - s0[2*i+1];
        vector<ZZ_p> ki(keybit.begin()+2*i, keybit.begin()+2*(i+1));
        ti = Fmult(ki, sminus);
        t.insert(t.end(), ti.begin(), ti.end());
    }
    for (int m=0; m<keybit.size(); m++)
    {
        p.push_back(s0[m]+t[m]);
    }
    return p;
}

vector<ZZ_p> mpc::inversePerm(vector<ZZ_p> pi)
{
    assertSize(pi, "inversePerm");
    vector<uint64_t> arr2(pi.size());
    vector<ZZ_p> inverse(pi.size());
    for (uint64_t i = 0; i < pi.size(); i++)
        arr2[conv<uint64_t>(pi[i]) - 1] = i + 1;
 
    for (uint64_t i = 0; i < pi.size(); i++)
    {
        inverse[i] = conv<ZZ_p>(arr2[i]);
    }
    return inverse;
}

void mpc::apply_perm_local(bool participate, vector<ZZ_p> &v, vector<ZZ_p> &pi)
{
    if (participate)
    {
        for (ZZ_p z : pi) 
        {
            uint64_t value = conv<uint64_t>(z);
            if (value > pi.size()) 
            {
                cout << "apply_perm_local." << endl;
                cout << "value: " << value << endl;
                cout << "size of vector: " << pi.size() << endl;
                throw logic_error("permutation is bigger than vector size.");
            }
        }
        if (pi.size() != v.size()/2)
            throw std::invalid_argument("local permutation should be half the size of shared vector");
        if (v.size() % 2 != 0)
            throw invalid_argument("v size is not divisible by 2.");
        vector<ZZ_p> v2(v.size());
        for (uint64_t i = 0; i < v.size()/2; i++)
        {
            if (2 * i + 1 >= v.size()) 
            {
                cout << "i too high" << i << "size: " << v.size() << endl;
            }
            try
            {
                uint64_t idx1 = conv<uint64_t>((pi[i])-1)*2;
                uint64_t idx2 = conv<uint64_t>((pi[i])-1)*2+1;
                v2[idx1] =v[2*i];
                v2[idx2]=v[2*i+1];
            }
            catch (const std::exception &e)
            {
                cerr << "Exception caught: " << e.what() << std::endl;
            }
        }
        // return v2;
        v = std::move(v2);
    }
}
void mpc::apply_perm_localM(bool participate, vector<vector<ZZ_p>> &v, vector<ZZ_p> &pi)
{
    if (participate)
    {
        for (ZZ_p z : pi) 
        {
            uint64_t value = conv<uint64_t>(z);
            if (value > pi.size()) 
            {
                cout << "apply_perm_local." << endl;
                cout << "value: " << value << endl;
                cout << "size of vector: " << pi.size() << endl;
                throw logic_error("permutation is bigger than vector size.");
            }
        }
        if (pi.size() != v.size())
            throw std::invalid_argument("local permutation should be half the size of shared vector");
        // if (v.size() % 2 != 0)
        //     throw invalid_argument("v size is not divisible by 2.");
        vector<vector<ZZ_p>> v2(v.size(), vector<ZZ_p>(v[0].size()));
        for (size_t i=0; i<v.size(); i++)
        {
            uint64_t idx = conv<uint64_t>(pi[i]-1);
            v2[idx] = v[i];
        }
        v = std::move(v2);
    }
}
void mpc::center_normalize_pheno() {
    double pheno_sum=0.0;
    double sq_error_sum = 0.0;
    if(this->pid == 0){
        double recv_sum, sq_error;
        this->dataowner.recv(recv_sum); //CP1 aggregates sum from each data owner
        pheno_sum += recv_sum;
        this->toOwner.send(pheno_sum);

        this->dataowner.recv(sq_error);
        sq_error_sum += sq_error; //CP1 aggregates sum of squared error
        this->toOwner.send(sq_error_sum);
    }
}
void mpc::permutPheno(int permut)
{

    vector<uint64_t> normalized_pheno;
    this->dataowner.recv(normalized_pheno);
    vector<ZZ_p> pheno = convVec(normalized_pheno);
    // ZZ_p mean = center_normalize_pheno(pheno);
    this->permutMat.push_back(pheno);
    for (size_t i=0; i<permut; i++)
    {
        vector<ZZ_p> permuted(pheno);
        vector<ZZ_p> pi=Frand(pheno.size()/2);
        shuffle(pi,permuted);
        this->permutMat.push_back(permuted);
    }
    if (this->pid == 0)
    {
        cout << "permutation matrix complete.\n";
    }
    // reveal_matrix(this->permutMat, this->permutMat, "permutpheno");
}
void mpc::shuffle(vector<ZZ_p> &pi, vector<ZZ_p> &a)
{
    assertSize(pi, "Shuffle");
    if (pi.size() != a.size())
        throw std::invalid_argument("Your shares and pi size are different");
    vector<ZZ_p> pi_m(pi.size()/2), pi_p(pi.size()/2);
    for (int i=0; i<pi.size()/2; i++)
    {
        pi_m[i] = pi[2*i];
        pi_p[i] = pi[2*i+1];
    }
    if (this->pid == 0)
    {   
        apply_perm_local(true, a, pi_m);
        reshare(a, 1);
        apply_perm_local(true, a, pi_p);
        reshare(a, 2);
        apply_perm_local(false, a, pi);
        reshare(a, 0);
    }
    else if (this->pid == 1)
    {
        apply_perm_local(false, a, pi);
        reshare(a, 1);
        apply_perm_local(true,a, pi_m);
        reshare(a, 2);
        apply_perm_local(true,a, pi_p);
        reshare(a, 0);
    }
    else
    {
        apply_perm_local(true, a, pi_p);
        reshare(a, 1);
        apply_perm_local(false, a, pi);
        reshare(a, 2);
        apply_perm_local(true, a, pi_m);
        reshare(a, 0);
    }
}
void mpc::unshuffle(vector<ZZ_p> &pi, vector<ZZ_p> &b)
{
    assertSize(pi, "unshuffle");
    vector<ZZ_p> pi_m, pi_p;
    for (size_t i=0; i<pi.size()/2; i++)
    {
        pi_m.push_back(pi[2*i]);
        pi_p.push_back(pi[2*i+1]);
    }
    vector<ZZ_p> inv_m(pi_m.size()), inv_p(pi_p.size());
    inv_m = inversePerm(pi_m);
    inv_p = inversePerm(pi_p); 
    if (this->pid == 0)
    {
        apply_perm_local(false,b,inv_p);
        reshare(b,0);
        apply_perm_local(true,b,inv_p);
        reshare(b,2);
        apply_perm_local(true,b,inv_m);
        reshare(b,1);
    }
    else if (this->pid == 1)
    {
        apply_perm_local(true,b,inv_p);
        reshare(b, 0);
        apply_perm_local(true,b,inv_m);
        reshare(b, 2);
        apply_perm_local(false,b,inv_p);
        reshare(b, 1);
    }
    else
    {
        apply_perm_local(true,b,inv_m);
        reshare(b, 0);
        apply_perm_local(false,b,inv_p);
        reshare(b, 2);
        apply_perm_local(true,b,inv_p);
        reshare(b, 1);
    }
}

void mpc::shuffleM(vector<ZZ_p> &pi, vector<vector<ZZ_p>> &a)
{
    assertSize(pi, "Shuffle");
    // if (pi.size() != a.size())
    //     throw std::invalid_argument("Your shares and pi size are different");
    vector<ZZ_p> pi_m(pi.size()/2), pi_p(pi.size()/2);
    for (int i=0; i<pi.size()/2; i++)
    {
        pi_m[i] = pi[2*i];
        pi_p[i] = pi[2*i+1];
    }
    if (this->pid == 0)
    {   
        apply_perm_localM(true, a, pi_m);
        reshareM(a, 1);
        apply_perm_localM(true, a, pi_p);
        reshareM(a, 2);
        apply_perm_localM(false, a, pi);
        reshareM(a, 0);
    }
    else if (this->pid == 1)
    {
        apply_perm_localM(false, a, pi);
        reshareM(a, 1);
        apply_perm_localM(true,a, pi_m);
        reshareM(a, 2);
        apply_perm_localM(true,a, pi_p);
        reshareM(a, 0);
    }
    else
    {
        apply_perm_localM(true, a, pi_p);
        reshareM(a, 1);
        apply_perm_localM(false, a, pi);
        reshareM(a, 2);
        apply_perm_localM(true, a, pi_m);
        reshareM(a, 0);
    }
}
void mpc::reveal_matrix(vector<vector<ZZ_p>>& geno, vector<vector<ZZ_p>>& pheno,string name)
{
    // vector<vector<double>> sliced_geno;
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sliced_geno;
    sliced_geno.resize(geno.size(), geno[0].size()/2);
    for (size_t i=0; i<geno.size(); i++)
    {
        vector<ZZ_p> row = reveal(geno[i], false); 
        vector<double> unscaledrow;
        for (size_t j=0; j<row.size(); j++)
        {
            int64_t unshifted;
            if (conv<uint64_t>(row[j]) > this->p/2)
                unshifted = conv<int64_t>(row[j]) - this->p;
            else
                unshifted = conv<int64_t>(row[j]);
            double r_nom = static_cast<double>(unshifted)/pow(10,6);
            unscaledrow.push_back(r_nom);
            // r2row.push_back(r_nom*r_nom);

        }
        // sliced_geno.push_back(unscaledrow);
        Eigen::Map<Eigen::VectorXd>(sliced_geno.row(i).data(), unscaledrow.size()) = 
            Eigen::VectorXd::Map(unscaledrow.data(), unscaledrow.size());
    }
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> pheno_mat;
    pheno_mat.resize(pheno.size(), pheno[0].size()/2);
    for (size_t i=0; i<pheno.size(); i++)
    {
        vector<ZZ_p> row = reveal(pheno[i], false); 
        vector<double> unscaledrow;
        for (size_t j=0; j<row.size(); j++)
        {
            int64_t unshifted;
            if (conv<uint64_t>(row[j]) > this->p/2)
                unshifted = conv<int64_t>(row[j]) - this->p;
            else
                unshifted = conv<int64_t>(row[j]);
            double r_nom = static_cast<double>(unshifted)/pow(10,6);
            unscaledrow.push_back(r_nom);
            // r2row.push_back(r_nom*r_nom);

        }
        // pheno_mat.push_back(unscaledrow);
        Eigen::Map<Eigen::VectorXd>(pheno_mat.row(i).data(), unscaledrow.size()) = 
            Eigen::VectorXd::Map(unscaledrow.data(), unscaledrow.size());
    }
    Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mult_res;
    mult_res.resize(sliced_geno.rows(), pheno_mat.rows());
    auto start = chrono::high_resolution_clock::now();
    mult_res = sliced_geno*pheno_mat.transpose();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInminutes = totalduration.count()/60.0;
    if (this->pid ==0)
    {
        cout << string("\nPlaintext Eigen matmult took " + to_string(totaldurationInminutes)+" minutes.") << endl;
        // writeEigenToTSV(sliced_geno, string(name + "_pQTL_geno"));
        // writeEigenToTSV(pheno_mat, string(name + "_pQTL_pheno"));
        // writeEigenToTSV(mult_res, string(name + "_pQTL_mmresult"));
    }
        
}
void mpc::calc_corr(Logger& cislogger, Logger& nominalLogger)
{
    if (this->permutMat.empty())
    {
        throw invalid_argument("permutation matrix has not been made.");
    }
    if (this->pid ==0)
        cout << "Matmult execution... " << flush;
    auto start = chrono::high_resolution_clock::now();
    // reveal_matrix(this->geno, this->permutMat, "plaintext");
    // reveal_matrix(this->permutMat, "permutMat");
    vector<vector<ZZ_p>> matmultresult = matmult(this->geno, this->permutMat);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = end - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/60.0;
    if (this->pid==0)     
        cout << totaldurationInminutes << " minutes" << endl;
    vector<vector<double>> matcorr, r2_nominal;
    // retrieving r_nominal and r2_nominal
    for (size_t i=0; i<matmultresult.size(); i++)
    {
        vector<ZZ_p> row = reveal(matmultresult[i], false); 
        vector<double> unscaledrow, r2row;
        for (size_t j=0; j<row.size(); j++)
        {
            int64_t unshifted;
            if (conv<uint64_t>(row[j]) > this->p/2)
                unshifted = conv<int64_t>(row[j]) - this->p;
            else
                unshifted = conv<int64_t>(row[j]);
            double r_nom = static_cast<double>(unshifted)/pow(10,12);
            unscaledrow.push_back(r_nom);
            r2row.push_back(r_nom*r_nom);

        }
        matcorr.push_back(unscaledrow);
        r2_nominal.push_back(r2row);
    }
    //calculating statistics
    int dof = this->zscores.size()-2;
    
    if (this->pid ==0)
    {
        // cout << "WHOP"<< endl;
        vector<double> std_ratio;
        // writematrixToTSV(matpheno, int startrow, int endrow, const string& name)
        // writematrixToTSV(matpheno, "bloodpheno");
        // writematrixToTSV(matgeno, "bloodgeno");
        // writematrixToTSV(matcorr, "bloodcorr");
        this->dataowner.recv(std_ratio);
        // writeVectorToTSV(std_ratio, "blood_stdratio");
        int buffersize;
        this->dataowner.recv(buffersize);
        char buffer[buffersize];
        this->dataowner.recv(buffer, buffersize);
        string serializedVar(buffer, sizeof(buffer) - 1);
        vector<string> cisVar;
        istringstream iss(serializedVar);
        string token;
        while (getline(iss, token, ';')) {
            cisVar.push_back(token);
        }
        // cout <<string("received number of cis variants: "+ to_string(cisVar.size()-1)+"\n");
        // cout << string("\nreceived std_ratio first: "+to_string(std_ratio[0])+ ", 752: "+to_string(std_ratio[752])+"\n");
        // cout <<string("received std ratio size: "+to_string(std_ratio.size())+"\n");
        // cout << string("Matmult shape: "+to_string(matcorr.size())+","+to_string(matcorr[0].size())+"\n");
        // cout << "WHOP2"<< endl;
        //Writing nominal pass
        for (int v=0; v<cisVar.size()-1; v++)
        {
            string message = string(cisVar[0]+"\t"+cisVar[v+1]+"\t"+to_string(v)+"\t");
            // logger.log(string(cisVar[0]+"\t"+cisVar[v]+"\t"+to_string(v)+"\t"));
            double dof = this->permutMat[0].size()/2 - 2;
            double slope = matcorr[v][0]*std_ratio[v];
            double tstat2 = dof * r2_nominal[v][0] / (1- r2_nominal[v][0]);
            double slope_se = abs(slope)/sqrt(tstat2);
            double r2_nom = r2_nominal[v][0];
            double r_nom = matcorr[v][0];
            double pval = getPvalueFromTstat2(tstat2, dof);
            message += string(to_string(dof)+"\t"+to_string(r_nom)+"\t"+to_string(r2_nom)+"\t"+to_string(sqrt(tstat2))+
            "\t"+to_string(pval)+"\t"+to_string(slope)+"\t"+to_string(slope_se));
            nominalLogger.log(message);
        }
        // cout << "WHOP3"<< endl;
        vector<vector<double>> transposed(matcorr[0].size(), vector<double>(matcorr.size()));
        for (size_t i = 0; i < matcorr.size(); ++i) {
            for (size_t j = 0; j < matcorr[i].size(); ++j) {
                transposed[j][i] = matcorr[i][j];
            }
        }
        cout << string("transposed R2 matrix shape: " + to_string(transposed.size()) + ", "+to_string(transposed[0].size())+"\n");
        // Getting r2_perm (maximum corr for each permutation)
        vector<double> r_perm;
        vector<size_t> idx; // holds the max values for corr and all permutations
        int counter =0;
        for (const auto& column : transposed) {
           
            auto max_it = std::max_element(column.begin(), column.end(),[](double a, double b) {
                return abs(a) < abs(b);
            });
            size_t max_index = distance(column.begin(), max_it);
            if(counter!=0)
                r_perm.push_back(*max_it);
            idx.push_back(max_index);
            counter++;
        }
        //Compute mean and variance of p-values
        vector < double > permPvalues;
        double init_dof = this->permutMat[0].size()/2 - 2;
        double mean = 0.0, variance = 0.0, beta_shape1 = 1.0, beta_shape2 = 1.0;
        
        // cout << "init dof: "<<init_dof << endl;
        if (doublevariance(r_perm, doublemean(r_perm)) != 0.0) 
        {
            // cout << "learning degree of freedom" << endl;
            learnDF(r_perm, init_dof);
            //LOG.println("  * Effective degree of freedom = " + sutils::double2str(true_df, 4));
		}
        for (int c = 0 ; c < r_perm.size() ; c ++) permPvalues.push_back(getPvalue(r_perm[c], init_dof)); //dof should be changed to true_dof
        for (int pv = 0 ; pv < permPvalues.size() ; pv++) mean += permPvalues[pv];
        mean /= permPvalues.size();
        for (int pv = 0 ; pv < permPvalues.size() ; pv++) variance += (permPvalues[pv] - mean) * (permPvalues[pv] - mean);
        variance /= (permPvalues.size() - 1);
        //Estimate shape1 & shape2
        if (cisVar.size()-1 > 1 && mean != 1.0) {
            beta_shape1 = mean * (mean * (1 - mean ) / variance - 1);
            beta_shape2 = beta_shape1 * (1 / mean - 1);
            // cout << "mleBeta" << endl;
            if (cisVar.size()-1 > 10) mleBeta(permPvalues, beta_shape1, beta_shape2);	//ML estimate if more than 10 variant in cis
        }
        // print_vector(r2_perm);
        int ix = idx[0]; //index of highest variant for corr
        double dof = init_dof;
        double slope = matcorr[ix][0]*std_ratio[ix];
        double tstat2 = dof * r2_nominal[ix][0] / (1- r2_nominal[ix][0]);
        double slope_se = abs(slope)/sqrt(tstat2);
        double r2_value = r2_nominal[ix][0];
        double pval_nom = getPvalueFromTstat2(tstat2,this->permutMat[0].size()/2 - 2);
        double pval_true_df = getPvalueFromTstat2(tstat2, init_dof);
        size_t count = count_if(r_perm.begin(), r_perm.end(), [r2_value](double value) {
            return value*value >= r2_value;
        });
        double pval_beta = pbeta(pval_true_df,beta_shape1, beta_shape2,1,0);
        double pval_adj = static_cast<double>(count + 1) / (r_perm.size());
        string cis_message = string(cisVar[0]+"\t"+cisVar[ix+1]+"\t"+to_string(ix)+"\t");
        cis_message+= string(to_string(beta_shape1)+"\t"+to_string(beta_shape2)+"\t"+to_string(init_dof)+"\t"+to_string(pval_true_df)+"\t"+
        to_string(matcorr[ix][0])+"\t"+to_string(r2_value)+"\t"+to_string(sqrt(tstat2))+"\t"+to_string(pval_nom)+"\t"+
        to_string(slope)+"\t"+to_string(slope_se)+"\t"+to_string(pval_adj)+"\t"+to_string(pval_beta));
        cislogger.log(cis_message);
        cout << string("\n----------------------------\nGene: "+cisVar[0]+"\nVariant: "+cisVar[ix]+
        "\nidx: "+to_string(ix)+
        "\nstd_ratio: "+to_string(std_ratio[ix])+
        "\nr_nom: "+to_string(matcorr[ix][0]) +
        "\nr2_nom: "+to_string(r2_value)+
        "\ndof: "+to_string(dof)+"\nslope: "+to_string(slope)+"\ntstat: "+to_string(sqrt(tstat2))+"\nslope_se: "+to_string(slope_se)+
        "\npval_true: "+to_string(pval_true_df)+"\npval_nom: "+to_string(pval_nom)+
        "\npval_perm: "+to_string(pval_adj)+"\nbeta_shape1: "+to_string(beta_shape1)+"\nbeta_shape2: "+to_string(beta_shape2)+"\npval_beta: "+to_string(pval_beta)+
        +"\n----------------------------\n");
        
        // vector<double> agg_stat;
        // agg_stat.push_back(static_cast<double>(ix));
        // agg_stat.push_back(slope);
        // agg_stat.push_back(sqrt(tstat2));
        // agg_stat.push_back(slope_se);
        // agg_stat.push_back(pval_adj);
        int complete = 1;
        this->toOwner.send(complete);
        // cout <<string("everything so far complete.\n");
        
    }
    
    // vector<ZZ_p> revealp = reveal(this->pheno);
    // if(this->pid !=0)
    //     print_vector(this->pheno);
    
}
void mpc::apply_shared_perm(vector<ZZ_p> &rho, vector<ZZ_p> &k)
{
    if (rho.size() != k.size()){
        cout << string("rho: "+ to_string(rho.size()) + ", k: "+ to_string(k.size()));
        throw std::invalid_argument("rho and k size don't match");
    }
        
    vector<ZZ_p> pi=Frand(k.size()/2);
    
    shuffle(pi, rho);
    shuffle(pi, k);
    vector<ZZ_p> reconstructed = reveal(rho, false); 
    assertSize(reconstructed, "apply_shared_perm");
    if (reconstructed.size() == k.size()/2)
    {
        vector<ZZ_p> v2(k.size());
        for (size_t i = 0; i < k.size()/2; i++)
        {
            uint64_t pi_i = conv<uint64_t>(reconstructed[i]);
            v2[(pi_i-1)*2] =k[2*i];
            v2[(pi_i-1)*2+1]=k[2*i+1];
        }
        k = std::move(v2);

    }
    else 
    {
        cout << "reconstructed size: " << reconstructed.size() << endl;
        cout << "vector size:" << k.size() << endl;
        throw logic_error("Apply shared perm: Asking to apply different size permutations");
    }
    // apply_perm_local(k, reconstructed);
    // print_vector(k);
}

void mpc::compose(vector<ZZ_p> &sigma, vector<ZZ_p> &rho)
{
    vector<ZZ_p> pi=Frand(sigma.size()/2);
    shuffle(pi,sigma);
    vector<ZZ_p> reconstructed = reveal(sigma, false);
    assertSize(reconstructed, "compose");
    vector<ZZ_p> sigma_inv = inversePerm(reconstructed);

    if (sigma_inv.size() == rho.size()/2)
    {
        vector<ZZ_p> v2(rho.size());
        for (size_t i = 0; i < rho.size()/2; i++)
        {
            uint64_t pi_i = conv<uint64_t>(sigma_inv[i]);
            v2[(pi_i-1)*2] =rho[2*i];
            v2[(pi_i-1)*2+1]=rho[2*i+1];
        }
        rho = std::move(v2);
    }
    else 
    {
        throw logic_error("Compose: Asking to apply different size permutations");
    }
    unshuffle(pi, rho);
    // return rho;
}
template <typename T>
void writeVectorToCSV(const std::vector<T>& data, int pid, int row, string name)
{
    std::string filename = "/gpfs/commons/groups/gursoy_lab/aychoi/eqtl/mpc/securesort/output/" + name + "_pid" + std::to_string(pid) + "_row" + to_string(row)+ ".csv";
    std::ofstream file(filename);
    if (file.is_open())
    {
        for (const T& value : data)
        {
            file << value << ",";
        }
        file.close();
        std::cout << "Vector successfully written to CSV file." << std::endl;
    }
    else
    {
        std::cout << "Error opening the file." << std::endl;
    }
}

void mpc::genperm(int row, string norm_method)
{
    if (this->pid ==0)
        cout << "Secure sort execution... ";
    auto sort_start = chrono::high_resolution_clock::now();
    vector<ZZ_p> k_i(2*this->n);
    k_i.assign(this->shares.begin(), this->shares.begin() + 2*this->n);
    // if (this->pid == 0)
    //     cout <<"assigned" <<endl;
    vector<ZZ_p> inv_rho, inv_cdf;
    inv_rho.assign(this->identity.begin(), this->identity.end());
    inv_cdf.assign(this->zscores.begin(), this->zscores.end());
    try
    {
        vector<ZZ_p> sigma = genbitperm(k_i);
        vector<ZZ_p> rho(sigma.size());
        for (int i=1; i<this->lk; i++)
        {
            k_i.assign(this->shares.begin() + 2*this->n * i, this->shares.begin() + 2*this->n * (i + 1));
            vector<ZZ_p> prev_sigma(sigma);
            apply_shared_perm(sigma, k_i);
            rho = genbitperm(k_i);
            compose(prev_sigma, rho);
            swap(sigma,rho);    
        }
        // if (this->pid == 0)
        //     cout << "loop finished.\n";
        vector<ZZ_p> reconstructed = reveal(sigma, false); 
        vector<uint64_t> rank_result = convVec(reconstructed);
        // if (this->pid ==0)
        //     cout << "\trank: " << string(to_string(rank_result[0]) +"\t" + to_string(rank_result[1]) + "\t" +to_string(rank_result[2]) +"\t" +to_string(rank_result[3])) << endl;
        apply_shared_perm(sigma, inv_rho);
        apply_shared_perm(inv_rho, inv_cdf);
        auto sort_end = chrono::high_resolution_clock::now();
        chrono::duration<double> sortduration = sort_end - sort_start;
        double totaldurationInSeconds = sortduration.count();
        double totaldurationInminutes = totaldurationInSeconds/60.0;
        if (this->pid == 0)
            cout << string(to_string(totaldurationInminutes) +" minutes.\n");
        // this->permutMat.push_back(inv_cdf);
        // for (size_t i=0; i<permut; i++)
        // {
        //     vector<ZZ_p> permuted(inv_cdf);
        //     vector<ZZ_p> pi=Frand(inv_cdf.size()/2);
        //     shuffle(pi,permuted);
        //     this->permutMat.push_back(permuted);
        // }
    }
    catch (const std::exception &e)
    {
        cerr << "Exception caught: " << e.what() << std::endl;
    }
    vector<uint64_t> invcdf_result = convVec(inv_cdf);
    // int complete = 1;
    this->toOwner.send(invcdf_result);
    
}
void mpc::receivePheno()
{
    vector<uint64_t> mat2;
    this->dataowner.recv(this->shape);
    this->dataowner.recv(mat2);
    for (int j=0; j<this->shape[0]; j++)
    {
        vector<ZZ_p> row;
        for (int k=0; k<this->shape[1]; k++)
        {
            row.push_back(to_ZZ_p(mat2[2*j*this->shape[1]+2*k]));
            row.push_back(to_ZZ_p(mat2[2*j*this->shape[1]+2*k+1]));
        }
        this->pheno.push_back(row);
    }
}
vector<vector<ZZ_p>> transpose(vector<vector<ZZ_p>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size()/2;
    vector<vector<ZZ_p>> transposed(cols, vector<ZZ_p>(rows*2));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transposed[j][2*i] = matrix[i][2*j];
            transposed[j][2*i+1] = matrix[i][2*j+1];
        }
    }
    return transposed;
}
void mpc::logRatio()
{
    // pseudo-reference is average log counts across all samples, for each gene
    uint64_t gene = this->pheno.size();
    uint64_t sample = this->pheno[0].size();
    // cout << string(to_string(this->pheno.size())+" / "+ to_string(this->pheno[0].size())+"\n");
    vector<uint64_t> sumvec;
    vector<double> pseudoref;
    vector<vector<double>> ratios(this->pheno.size(), vector<double>(this->pheno[0].size()));
    vector<double> sendback;
    for (int i=0; i<gene; i++)
    {
        ZZ_p sum=to_ZZ_p(0);
        for (int j=0; j<sample; j++)
        {
            sum+=this->pheno[i][j];
        }

        // double avg = conv<uint64_t>(sum)/static_cast<double>(sample);
        // cout <<sum<<endl;
        sumvec.push_back(conv<uint64_t>(sum));
    }
    // print_vector(sumvec);
    vector<uint64_t> v1, v2;
    this->toMinus.send(sumvec);
    this->toPlus.send(sumvec);
    this->fromMinus.recv(v1);
    this->fromPlus.recv(v2);
    for (int i=0; i<sumvec.size();i++)
    {
        ZZ_p totalsum=to_ZZ_p(sumvec[i]);
        totalsum+=to_ZZ_p(v1[i]);
        totalsum+=to_ZZ_p(v2[i]);
        // cout <<totalsum;
        pseudoref.push_back(conv<uint64_t>(totalsum)/static_cast<double>(sample));
    }
    // vector<ZZ_p> pseudosum = reveal(pseudoref,false);
    // print_vector(pseudoref);
    // for (int i=0; i<gene; i++)
    // {
    //     for(int j=0; j<sample; j++)
    //     {
    //         ratios[i][j] = conv<uint64_t>(this->pheno[i][j]) - pseudoref[i];
    //         sendback.push_back(ratios[i][j]);
    //     }
    // }
    // print_vector(sendback);
    this->toOwner.send(pseudoref);
}
void mpc::center_normalize_geno()
{
    if(this->pid ==0)
    {
        vector<double> data_sum;
        this->dataowner.recv(data_sum);
        vector<double> geno_sum(data_sum.size(),0.0);
        for (int i=0; i<data_sum.size(); ++i)
        {
            geno_sum[i] += data_sum[i];
        }
        this->toOwner.send(geno_sum);
        vector<double> error_sum;
        this->dataowner.recv(error_sum);
        vector<double> geno_error_sum(error_sum.size(),0.0);
        for (int i=0; i<error_sum.size(); ++i)
        {
            geno_error_sum[i] += error_sum[i];
        }
        this->toOwner.send(geno_error_sum);
    }

}
void mpc::receiveGeno()
{
    vector<uint64_t> mat1, genoshape;
    this->dataowner.recv(genoshape);
    // this->shape.push_back(genoshape);
    this->dataowner.recv(mat1);
    for (int j=0; j<genoshape[0]; j++)
    {
        vector<ZZ_p> row;
        for (int k=0; k<genoshape[1]; k++)
        {
            row.push_back(to_ZZ_p(mat1[2*j*genoshape[1]+2*k]));
            row.push_back(to_ZZ_p(mat1[2*j*genoshape[1]+2*k+1]));
        }
        this->geno.push_back(row);
    }
}
vector<vector<ZZ_p>> mpc::matmult(vector<vector<ZZ_p>>& mat1, vector<vector<ZZ_p>>& mat2)
{
    
    if (mat1[0].size() != mat2[0].size()) 
        throw logic_error("mat1 column size and mat2 row size do not match.");
    // vector<vector<ZZ_p>> result(mat1.size(), vector<ZZ_p>(mat2.size()));
    // vector<vector<ZZ_p>> transposed = transpose(mat2);
    auto minus_it = this->seedpair.find((this->pid + 2) % 3);
    auto plus_it = this->seedpair.find((this->pid + 1) % 3);
    PRNG &minus_PRNG = *(minus_it->second);
    PRNG &plus_PRNG = *(plus_it->second);
    // vector<uint64_t> to_send;
    const int nnRows = mat1.size();
    const int nnCols = mat1[0].size() / 2;

    MatrixXi mat1_share1, mat1_share2, mat2_share1, mat2_share2;
    mat1_share1.resize(mat1.size(), mat1[0].size()/2);
    mat1_share2.resize(mat1.size(), mat1[0].size()/2);
    mat2_share1.resize(mat2[0].size()/2, mat2.size());
    mat2_share2.resize(mat2[0].size()/2, mat2.size()); // split into transposes
    ZZtoEigen(mat1, mat1_share1,mat1_share2,false);
    ZZtoEigen(mat2, mat2_share1,mat2_share2,true);
    
    MatrixXi result1;
    result1.resize(mat1_share1.rows(), mat2_share1.cols());
    // MatrixXd result1(mat1_share1.rows(), mat2_share1.cols(), Eigen::RowMajor);
    if (this->pid ==0)
    {
        cout << "Eigen matmult..." << flush;
    }
        
    auto start = chrono::high_resolution_clock::now();
    result1.noalias() = mat1_share1*mat2_share1 + mat1_share1*mat2_share2 + mat1_share2*mat2_share1;
    auto matmult = chrono::high_resolution_clock::now();
    chrono::duration<double> matmult_dur = matmult - start;
    double matmult_InSeconds = matmult_dur.count();
    double matmult_Inminutes = matmult_InSeconds/ 60.0;
    if (this->pid ==0)
        cout << string("finished "+to_string(matmult_Inminutes)+" minutes. Creating random number matrices...") << flush;
    MatrixXi random1, random2;
    random1.resize(result1.rows(), result1.cols());
    random2.resize(result1.rows(), result1.cols());
    random1.setConstant(minus_PRNG.get<uint64_t>()%this->p);
    random2.setConstant(plus_PRNG.get<uint64_t>()%this->p);
    if (this->pid ==0)
        cout << "finished. Adding randomness..." << flush;
    result1 = result1+random1-random2;

    if (this->pid ==0)
        cout << "finished. Flattening to send vector..." << flush;
    vector<uint64_t> to_send(result1.data(), result1.data()+result1.size());

    if (this->pid ==0)
        cout << string("finished with vector size:"+ to_string(to_send.size()) +".\nSending and receiving...") << endl;
    
    //////////
    const size_t chunkSize = 200000;
    vector<uint64_t> received_data;
    for (size_t i = 0; i < to_send.size(); i += chunkSize) {
        size_t remaining = min(chunkSize, to_send.size() - i);
        vector<uint64_t> chunk_send(to_send.begin() + i, to_send.begin() + i + remaining);
        // Send 'to_send' in chunks of size 500 (or less if remaining is less than 500)
        this->toMinus.send(chunk_send);
        vector<uint64_t> to_receive;
        this->fromPlus.recv(to_receive);
        received_data.insert(received_data.end(), to_receive.begin(), to_receive.end());
        if(this->pid == 0)
            cout << "\tSent and received " << remaining << " elements in this chunk." << endl;
    }

    if (this->pid ==0)
        cout << "Send and receive done:" << flush;
    auto postRecv = chrono::high_resolution_clock::now();
    chrono::duration<double> totalduration = postRecv - start;
    double totaldurationInSeconds = totalduration.count();
    double totaldurationInminutes = totaldurationInSeconds/ 60.0;
    if (this->pid==0){
        cout << totaldurationInminutes << " minutes" << endl;
        // cout << string(to_string(to_send.size())+", "+to_string(received_data.size()))<<endl;
    }     

    vector<vector<ZZ_p>> finalresult(mat1.size(), vector<ZZ_p>(mat2.size()*2));
    EigentoZZ(to_send, received_data, finalresult);

    return finalresult;
}
int mpc::validgene()
{
    int returned;
    this->dataowner.recv(returned);
    return returned;
}
void mpc::close()
{
    // std::lock_guard<std::mutex> lock(mtx);
    this->dataowner.close();
    this->fromPlus.close();
    this->fromMinus.close();
    this->toPlus.close();
    this->toMinus.close();
    if (this->pid == 0)
        cout << "channels closed.\n";
}