#ifndef BLS_ZKPoS_H
#define BLS_ZKPoS_H
#define BITSIZE 160
#define BLOCKSIZE (8192*8192*4)
#define SAMPLESIZE 1000
#include <pbc.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <set>
#include <chrono>
#include <iterator>
#include "common.h"

class Proof
{
public:
    element_t u, lambda, ibzl;
    Proof(pairing_t pairing) {element_init_Zr(u, pairing); element_init_G1(lambda, pairing); element_init_Zr(ibzl, pairing);}
    ~Proof() {}
};

class safe_mpz
{// mpz_t has not a constructor, cannot use vector
public:
    mpz_t z;
    safe_mpz(){mpz_init(z);}
    ~safe_mpz(){//mpz_clear(z);
    }
};

class myelement
{
public:
    element_t e;
    myelement(){};
    ~myelement(){};
};

class BLS_ZKPoS {
public:
    //int keyGen(mpz_t security_param);
    pairing_t pairing;

    BLS_ZKPoS()
    {
        char param[1024];
        FILE* pParam=fopen("../param/a.param", "r");
        size_t count = fread(param, 1, 1024, pParam);
        fclose(pParam);
        if (!count) pbc_die("no param");
        pairing_init_set_buf(this->pairing, param, count);
        element_init_G1(this->u, this->pairing);
        element_init_G1(this->g1, this->pairing);
        element_init_G1(this->hash_value, this->pairing);
        element_init_Zr(this->x, this->pairing);
        element_init_Zr(this->rm, this->pairing);
        element_init_Zr(this->rsigma, this->pairing);
        element_init_Zr(this->p, this->pairing);
        element_init_Zr(this->commitment, this->pairing);
        element_init_Zr(this->small_hash_value, this->pairing);
        element_init_G2(this->g, this->pairing);
        element_init_G2(this->v, this->pairing);
        element_init_GT(this->pke, this->pairing);
        element_init_GT(this->R, this->pairing);

    }
    ~BLS_ZKPoS()
    {

    }
    int keyGen();// security param 1^k
    int sigGen(std::vector<safe_mpz> file, std::vector<myelement>& auth, std::vector<myelement>& names);
//    int sig(const element_t sk, const element_t message, element_t& sig);
//    int sigVerify(const element_t sk, const element_t sig);
    int prove(std::vector<int> index, std::vector<myelement> v, std::vector<safe_mpz> files, std::vector<myelement> auth, Proof& pi);
    int commit();
    int challenge(std::vector<int>& index, std::vector<myelement>& v, int len);
    int verify(std::vector<int> index, Proof pi, std::vector<myelement> names, std::vector<myelement> v);

    int exportPk(std::string pkFileName = std::string("rsa_test.pk"));
    int exportSk(std::string skFileName = std::string("rsa_test.sk"));
    int exportKeys(); // both pk and sk to two files


private:
    int H(mpz_t z);//G1
    int h(element_t R);
    element_t ssk, spk, x, v, g, g1, u, pke, hash_value, rm, rsigma, p, commitment, R, small_hash_value;
};


#endif //BLS_ZKPoS_H
