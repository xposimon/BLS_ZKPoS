cmake_minimum_required (VERSION 2.6)
project (BLS_ZKPoS)
include_directories(/usr/local/include/pbc/)
include_directories(/usr/local/lib/)

add_executable(bls_test main.cpp BLS_ZKPoS.cpp BLS_ZKPoS.h common.cpp)
target_link_libraries(bls_test gmp)
target_link_libraries(bls_test pbc)