#include <cstdlib>
#include "../../LaiCa/LaiCa.h"
#include <iomanip>
#include <itpp/comm/modulator.h>
#include<iostream>
#include<vector> //Vectors
#include<stdlib.h>
#include<string>
#include <fstream> //I/O operations
#include <math.h> //For pow()
#include<time.h>
#include<stdint.h>
#include<complex.h>

using namespace std;
using namespace laica;
using namespace itpp;

int main(int argc, char** argv) 
{
    //Open results file
    ofstream outputfile("results.txt");
    // If we couldn't open the output file stream for writing 
    if (!outputfile)
    {
            // Print an error and exit
            cout << "results.txt could not be opened for writing!" << endl;
            exit(1);
    }
    
    int M = 4; // Modulation cardinality 
    int bitSym = log2(M);
    
    // Code
    int n_LDPC = 64800;
    const int num_r = 1; //code rate numerator
    const int den_r = 2; //code rate denominator 
    const int maxItLDPC = 10; // iterations of LDPC decoder
    const int maxNumErr = 100;
    const float R = 1.0 / 2.0;
    const int nSymb = n_LDPC / log2(M);
    
    // Code
    DVBS2X_Code ldpc;
    const double rateLDPC = ldpc.set_parameters(n_LDPC, num_r, den_r);
    const int k_LDPC = R * n_LDPC;
    
    // PSK    
    laica::PSK psk;
    psk.set_parameters(M);    
    
    // Detector
    SBS_Detector sbs;
    sbs.set_parameters(psk.get_symbols(), 1.0);
    
    // Mapper
    Mapper map;
    map.set_parameters(psk.get_bits2symbols(), M, n_LDPC);
    
    //Ber meter
    BER_Meter meter;
    meter.set_parameters(k_LDPC);
    
    vec llr_ext;
    vec soft_dec;
    
    vec EbN0 = "0:1";
    
    for (int l = 0; l < size(EbN0); l++)
    {
        const double noisevar = 1 / (R * inv_dB(EbN0(l)) * log2(M));
        const double phasevar = 1 / inv_dB(0);
        meter.reset();
    
//        while (meter.get_frame_errors() < maxNumErr) 
//        {
            bvec msg = randb(k_LDPC);
            bvec c = ldpc.encode(msg);
            
            cvec symb(n_LDPC / bitSym);
            symb = psk.modulate_bits(c);    
            
            cvec x = symb + randn_c(nSymb) * sqrt(noisevar);
            
            // Phase noise distribution
            vec deltaTheta = randn(nSymb) * sqrt(phasevar);
            vec theta(nSymb);
            theta(0) = randu() * 2 * pi;
            for (int i = i; i < nSymb; i++)
            {
                theta(i) = theta(i-1) + deltaTheta(i);
            }
            
            cvec r(nSymb);
            for(int i = 0; i < nSymb; i++)
            {
                r(i) = symb(i) * exp(complex<double>(0,theta(i))) + randn_c() * sqrt(noisevar);
                outputfile << i << " " << r(i) << endl;
            }
            
            
            
//            mat pr_soft;
//            sbs.soft_detection(x, pr_soft);
//            
//            map.reset();
//            vec llr = map.softsymbols_to_llr(pr_soft);
            
            //        ldpc.reset();
            //        int wsyn = 1;
            //        for (int it = 0; it < maxItLDPC; it++)
            //        {
            //            ldpc.decode(llr, soft_dec, llr_ext);
            //            wsyn = ldpc.weight_syndrome(soft_dec);
            //            if (wsyn == 0)
            //                break;
            //        }
            //        if (wsyn == 0)
            //            break;
            //        
            //        meter.run(msg, soft_dec.get(0, k_LDPC - 1));
            //        cout << meter.get_BER() << " " << meter.get_frame_errors() << endl;
//        }    
    }
    outputfile.close();
    
    return 0;
}

