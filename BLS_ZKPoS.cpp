#include "BLS_ZKPoS.h"

int BLS_ZKPoS::H(mpz_t z)
{
    // TODO simple Hash, not secure, return u^z
    element_mul_mpz(this->hash_value, this->u, z); // TODO not clear tmp, is there a destructor in mpz_t?
    return 1;
};

int BLS_ZKPoS::h(element_t R)
{
    unsigned char hash[1024];
    element_to_bytes(hash, R);

    element_from_hash(this->small_hash_value, hash, 1024);
}

int BLS_ZKPoS::exportPk(std::string pkFileName)
{
    FILE *pFile;
    pFile = fopen(pkFileName.c_str(), "w");
    return 1;
}

int BLS_ZKPoS::exportSk(std::string skFileName)
{
    FILE *pFile;
    pFile = fopen(skFileName.c_str(), "w");

    return 1;
}

int BLS_ZKPoS::exportKeys()
{
    this->exportPk();
    this->exportSk();
}

int BLS_ZKPoS::keyGen()
{
    element_random(this->g);
    element_random(this->g1);
    element_random(this->x);
    element_random(this->u);
    element_pow_zn(this->v, this->g, this->x);
    element_pairing(this->pke, this->u, this->v);
}

int BLS_ZKPoS::sigGen(std::vector<safe_mpz> file, std::vector<myelement>& auth, std::vector<myelement>& names)
{
    int len = file.size();
    auth.resize(len);
    names.resize(len);

    mpz_t mptmp;
    mpz_init(mptmp);
    std::string t1, t2;
    for (int i = 0; i < len; i++)
    {
        element_init_Zr(names[i].e, this->pairing);
        element_init_G1(auth[i].e, this->pairing);
        element_random(names[i].e);


        mpz_set_ui(mptmp, i);
        mp2bitString(mptmp, t2);

        element_to_mpz(mptmp, names[i].e);
        mp2bitString(mptmp, t1);




        t1 = t2+t1; // name||i bit concat
        std::reverse(t1.begin(), t1.end());

        mpz_set_str(mptmp, t1.c_str(), 2);

        gmp_printf("%Zd\n", mptmp);
        this->H(mptmp);

        element_pow_mpz(auth[i].e, this->u, file[i].z);
        element_mul(auth[i].e, this->hash_value, auth[i].e);
        element_pow_zn(auth[i].e, auth[i].e, this->x);
    }
    mpz_clear(mptmp);
    return 1;
}

int BLS_ZKPoS::commit()
{
    element_t tmp;
    element_init_GT(tmp, this->pairing);

    element_random(this->rm);
    element_random(this->rsigma);
    element_random(this->p);

    element_printf("rm:%B\n", this->rm);

    element_pairing(tmp, this->g1, this->g);
    element_pow_zn(tmp, tmp, this->rsigma);
    element_set(this->R, tmp);
    element_pairing(tmp, this->u, this->v);
    element_pow_zn(tmp, tmp, this->rm);
    element_mul(this->R, this->R, tmp);
    this->h(this->R);
    element_set(this->commitment, this->small_hash_value);
    element_printf("commit:%B\n", this->commitment);
    return 1;
}

int BLS_ZKPoS::challenge(std::vector<int>& index, std::vector<myelement>& v, int len)
{
    // simulate challenge
    int samplesize = min(SAMPLESIZE, len);
    v.resize(samplesize);
    index.resize(samplesize);
    std::default_random_engine random;
    std::uniform_int_distribution<int> u(0, len-1);
    std::set<int> chosen;
    int tmpc;
    for( int i = 0; i < samplesize; i++)
    {
        element_init_Zr(v[i].e, this->pairing);
        element_random(v[i].e);
        if (samplesize < SAMPLESIZE)
            index[i] = i;

        //mpz_set_ui(coeff[i].z, 1);
    }
    if (samplesize < SAMPLESIZE)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle (index.begin(), index.end(), std::default_random_engine(seed));
    }
    else {
        for (int i = 0; i < samplesize; i++) {
            tmpc = u(random);
            while (chosen.find(tmpc) != chosen.end())
                tmpc = u(random);
            index[i] = tmpc;
            chosen.insert(tmpc);
        }
    }

    //mpz_set_ui(R, 1);

    return 1;
}

int BLS_ZKPoS::prove(std::vector<int> index, std::vector<myelement> v, std::vector<safe_mpz> files, std::vector<myelement> auth, Proof& pi)
{
    int len = index.size();
    element_t up, sigma;
    element_init_Zr(up, this->pairing);
    element_init_G1(sigma, this->pairing);
    element_set0(up);
    element_set1(sigma);
    element_t tmp;
    element_init_Zr(tmp, this->pairing);
    for (int i =0;i<len; i++)
    {
        element_mul_mpz(tmp, v[i].e, files[index[i]].z);
        element_printf("tmp:%B %B\n", v[i].e, tmp);
        element_add(up, up, tmp); // u=sum(vi mi)
    }
    element_init_G1(tmp, this->pairing);
    for (int i = 0; i < len; i++)
    {
        element_pow_zn(tmp, auth[index[i]].e, v[i].e);
        element_mul(sigma, sigma, tmp);
        element_printf("sigma: %B\n %B\n", sigma, tmp);
    }
    element_mul(pi.u, this->commitment, up);
    element_add(pi.u, pi.u, this->rm);
    element_pow_zn(pi.lambda, this->g1, this->p);
    element_mul(pi.lambda, sigma, pi.lambda);
    element_mul(pi.ibzl, this->commitment, this->p);
    element_add(pi.ibzl, this->rsigma, pi.ibzl);

    element_t T1, T2;
    element_init_G1(T1, this->pairing);
    element_init_GT(T2, this->pairing);
    element_pow_zn(T1, this->g1, this->p);
    element_mul(T1, T1, sigma);
    element_pow_zn(T1, T1, this->commitment);
    element_pairing(T2, T1, this->g);

    return 1;
}

int BLS_ZKPoS::verify(std::vector<int> index, Proof pi, std::vector<myelement> names, std::vector<myelement> v)
{
    element_t left, right, tmp;
    element_t g1tmp, g1tmp2;
    element_init_G1(g1tmp2, this->pairing);
    element_init_G1(g1tmp, this->pairing);
    element_init_GT(left, this->pairing);
    element_init_GT(right, this->pairing);
    element_init_GT(tmp, this->pairing);
    element_set1(g1tmp2);

    element_pow_zn(g1tmp, pi.lambda, this->commitment);
    element_pairing(left, g1tmp, this->g);

    element_mul(left, this->R, left);

    element_t T1, T2;
    element_init_GT(T1, this->pairing);
    element_init_GT(T2, this->pairing);
    element_pairing(T1, this->g1, this->g);
    element_pow_zn(T1, T1, this->rsigma);
    element_pairing(T2, this->u, this->v);
    element_pow_zn(T2, T2, this->rm);
    element_mul(T1, T1, T2);

    element_pairing(right, this->g1, this->g);
    element_pow_zn(right, right, pi.ibzl);
    int len = index.size();
    mpz_t mptmp;
    mpz_init(mptmp);
    std::string t1, t2;

    for(int i = 0; i < len; i++)
    {
        mpz_set_ui(mptmp, index[i]);
        mp2bitString(mptmp, t2);

        element_to_mpz(mptmp, names[index[i]].e);
        mp2bitString(mptmp, t1);

        t1 = t2+t1; // name||i bit concat
        std::reverse(t1.begin(), t1.end());

        mpz_set_str(mptmp, t1.c_str(), 2);
        this->H(mptmp);
        element_pow_zn(g1tmp, this->hash_value, v[i].e);
        element_mul(g1tmp2, g1tmp2, g1tmp);
    }

    element_pow_zn(g1tmp2, g1tmp2, this->commitment);
    element_pow_zn(g1tmp, this->u, pi.u);
    element_mul(g1tmp, g1tmp, g1tmp2);
    element_pairing(tmp, g1tmp, this->v);
    element_mul(right, tmp, right);

    element_printf("res:%B\n %B\n", left, right);

    if(!element_cmp(left, right))
        return 1;
    return 0;

}