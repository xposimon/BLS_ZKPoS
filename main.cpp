#include <fstream>
#include "BLS_ZKPoS.h"
#define SECURITY_BITS 1024

int main()
{
    //std::cout<<"BLS_ZKPoS test codes"<<std::endl;
    BLS_ZKPoS* bls_test = new BLS_ZKPoS();
#ifdef TEST
    clock_t startTime, endTime;
    std::vector<double> stepTime;
    int file_length;
    startTime = clock();
#endif
    bls_test->keyGen();
#ifdef TEST
    endTime = clock();
    stepTime.push_back((double)(endTime-startTime)/CLOCKS_PER_SEC);
#endif
    bls_test->exportKeys();

    // bls_test->dbg_importPk(5, 5, 163, 2, 1024, 4761);

    std::cout<<"\nReading File: "<<std::endl;
    std::string filepath = std::string("../inputdata/1.in");
    std::ifstream fin(filepath.c_str(), std::ios::binary);

    std::vector<myelement> auth;
    std::vector<safe_mpz> fileBlocks;
    std::vector<myelement> names;
    std::vector<int> index;

    if (fin)
    {
        // file handle to vector<mpz_t> files with each mpz_t content of a block

        fin.seekg(0, fin.end);
        int length = (int)fin.tellg();
        file_length = length;
        std::cout<<"File length: "<<length<<std::endl;
        fin.seekg(0, fin.beg);
        std::string buffer(BLOCKSIZE, '\0');
        fileBlocks.resize(length%BLOCKSIZE==0?length/BLOCKSIZE:((int)(length/BLOCKSIZE)+1));

        int readSize = min(length,BLOCKSIZE);
        int count = 0, base = 1;

        while (readSize > 0 && fin.read(&buffer[0], readSize))
        {
            mpz_set_ui(fileBlocks[count].z, 0);
            base = 1;
            //std::cout<< (int)buffer[0]<< " "<<count<<std::endl;
            for (int i = 0; i < readSize; i++){
                mpz_add_ui(fileBlocks[count].z, fileBlocks[count].z, base*(int)buffer[i]);
                base *= 256;
            }

            length -= BLOCKSIZE;
            readSize = min(length, BLOCKSIZE);
            count ++;
        }

        auth.resize(count);
        names.resize(count);

        // generate tags
#ifdef TEST
        std::cout<<"Start authenticator generation"<<std::endl;
        startTime = clock();
#endif
        bls_test->sigGen(fileBlocks, auth, names);
#ifdef TEST
        endTime = clock();
        stepTime.push_back((double)(endTime-startTime)/CLOCKS_PER_SEC);
#endif
    }
    fin.close(); // finish file handle
    std::cout<<"commitment\n"<<std::endl;
    mpz_t commitment;
    mpz_init(commitment);
    bls_test->commit();

#ifdef DEBUG
    FILE *pFile=fopen("randoma.in", "w");
    fprintf(pFile, "a:");
    mpz_out_str(pFile, 10, commitment);
    fprintf(pFile, "\n");
#endif
    //simulate challenge
    std::vector<myelement> v;
    bls_test->challenge(index, v, fileBlocks.size());

#ifdef TEST
    startTime = clock();
#endif
    // generate proof
    Proof pi(bls_test->pairing);
    bls_test->prove(index, v, fileBlocks, auth, pi);
#ifdef TEST
    endTime = clock();
    stepTime.push_back((double)(endTime-startTime)/CLOCKS_PER_SEC);
#endif


#ifdef DEBUG
    fprintf(pFile, "t:");
    mpz_out_str(pFile, 10, pi.t);
    fprintf(pFile, "\n");
    fprintf(pFile, "sigma:");
    mpz_out_str(pFile, 10, pi.sigma);
    fprintf(pFile, "\n");
    fprintf(pFile, "u:");
    mpz_out_str(pFile, 10, pi.u);
    fprintf(pFile, "\n");
    fprintf(pFile, "R:");
    mpz_out_str(pFile, 10, R);
    fprintf(pFile, "\n");


    pFile=fopen("paraminfo.in", "w");
    fprintf(pFile, "file:");
    for (int i = 0; i < fileBlocks.size(); i++)
        mpz_out_str(pFile, 10, fileBlocks[i].z), fprintf(pFile, " ");
    fprintf(pFile,"\n");

    fprintf(pFile, "tags:");
    for (int i = 0; i < fileBlocks.size(); i++)
        mpz_out_str(pFile, 10, tags[i].z), fprintf(pFile, " ");
    fprintf(pFile,"\n");

    fprintf(pFile, "names:");
    for (int i = 0; i < fileBlocks.size(); i++)
        mpz_out_str(pFile, 10, names[i].z), fprintf(pFile, " ");
    fprintf(pFile,"\n");

    fprintf(pFile, "r:");
    for (int i = 0; i < fileBlocks.size(); i++)
        mpz_out_str(pFile, 10, r[i].z), fprintf(pFile, " ");
    fprintf(pFile,"\n");

    fprintf(pFile, "coeff:");
    for (int i = 0; i < fileBlocks.size(); i++)
        mpz_out_str(pFile, 10, coeff[i].z), fprintf(pFile, " ");
    fprintf(pFile,"\n");
#endif

#ifdef TEST
    startTime = clock();
#endif
    //simulate verify
    int res = bls_test->verify(index, pi, names, v);
#ifdef TEST
    endTime = clock();
    stepTime.push_back((double)(endTime-startTime)/CLOCKS_PER_SEC);
#endif
    std::cout<<"\nVerify result is "<<(bool)res<<std::endl;
#ifdef TEST
    std::cout<<"\nTime analysis:"<<std::endl;
    std::cout<<"File size:"<<file_length<<" bytes"<<std::endl;
    std::cout<<"Key generation time:"<<stepTime[0]<<" s"<<std::endl;
    std::cout<<"Tag generation time:"<<stepTime[1]<<" s"<<std::endl;
    std::cout<<"Proof generation time:"<<stepTime[2]<<" s"<<std::endl;
    std::cout<<"Verify time:"<<stepTime[3]<<" s"<<std::endl;
#endif
//    for(int i = 0; i < fileBlocks.size(); i++)
//        mpz_clear(fileBlocks[i]), mpz_clear(names[i]), mpz_clear(tags[i]), mpz_clear(coeff[i]);

    delete bls_test;
    return 0;
}