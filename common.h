//
// Created by nclab on 1/3/19.
//
//#define DEBUG
#define TEST

#ifndef BLS_ZKPOS_COMMON_H
#define BLS_ZKPOS_COMMON_H
#include <gmp.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
int min(int a, int b);

void mp2bits(const mpz_t z, std::vector<bool>& bits);

void mp2bitString(const mpz_t z, std::string& bitString);

#endif //BLS_ZKPOS_COMMON_H
