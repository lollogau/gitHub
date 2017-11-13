#include <cstdlib>
#include "../../LaiCa/LaiCa.h"
#include <iomanip>
#include <itpp/comm/modulator.h>

using namespace std;
using namespace laica;

int main(int argc, char** argv) 
{
    int M = 8; // Modulation cardinality 
    int bitSym = log2(M);

    // Code
    int n_LDPC = 64800;
    const int num_r = 1; //code rate numerator
    const int den_r = 2; //code rate denominator 
    const int maxItLDPC = 10; // iterations of LDPC decoder
    const int maxNumErr = 30;
    const float R = 1.0 / 2.0;
    const int nSymb = n_LDPC / log2(M);
    

    return 0;
}

