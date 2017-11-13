#include <cstdlib>
#include "../../LaiCa/LaiCa.h"
#include <iomanip>
#include <itpp/comm/modulator.h>

using namespace std;
using namespace laica;
using namespace itpp;

int main(int argc, char** argv) 
{
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
    
    for (float EbN0 = 0; EbN0 < 5; EbN0++) 
    {
        const double noisevar = 1 / (R * inv_dB(EbN0) * log2(M));
        meter.reset();
        
        while (meter.get_frame_errors() < maxNumErr) 
        {
            bvec msg = randb(k_LDPC);
            bvec c = ldpc.encode(msg);
            
            cvec symb(n_LDPC / bitSym);
            symb = psk.modulate_bits(c);    

            cvec x = symb + randn_c(nSymb) * sqrt(noisevar);
            
            mat pr_soft;
            sbs.soft_detection(x, pr_soft);
                        
            map.reset();
            vec llr = map.softsymbols_to_llr(pr_soft);
            
            ldpc.reset();
            int wsyn = 1;
            for (int it = 0; it < maxItLDPC; it++)
            {
                ldpc.decode(llr, soft_dec, llr_ext);
                wsyn = ldpc.weight_syndrome(soft_dec);
                if (wsyn == 0)
                    break;
            }
            if (wsyn == 0)
                    break;
            
            meter.run(msg, soft_dec.get(0, k_LDPC - 1));
            cout << meter.get_BER() << " " << meter.get_frame_errors() << endl;
        }
        cout << EbN0 << " " << meter.get_BER() << " ----------------------------- " << endl;
    }
    
    
    return 0;
}

