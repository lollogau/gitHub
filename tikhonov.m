clc
clear all

M = 4; % COnstellation Cardinality

hEnc = comm.LDPCEncoder;
hMod = comm.PSKModulator(M, 'BitInput',true, 'PhaseOffset', 0);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',5);
hDemod = comm.PSKDemodulator(4, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance', 1/10^(hChan.SNR/10));
hDec = comm.LDPCDecoder('OutputValue','Whole codeword','DecisionMethod','Soft decision');
hError = comm.ErrorRate;

constSymb = constellation(hMod); % Constellation
constBits = zeros(M, log2(M)); % Vector containing the bits of constellation symbols
for i = 0:2^(log2(M))-1
    symbol = hMod(double(dec2bin(i, log2(M))).' - 48);
    index = find(constSymb == symbol);
    constBits(index,:) = double(dec2bin(i, log2(M))) - 48;
    constSymbMod(i+1) = symbol;
end
release(hMod);

noisevar = 1 / 10^(hChan.SNR/10);
phasevar = 6;

data = logical(randi([0 1], 32400, 1)); % Message bits
encodedData = hEnc(data); % Encoded bits
c = hMod(encodedData); % Modulated symbols
nSymb = length(c); % Number of Symbols

% Pilot symbols insertion
pilotStep = 10;
nPilot = length(c) / pilotStep ; % Number of pilot symbols
cPilot = zeros(1, nSymb +nPilot);
j = 1;
for i = 1:nSymb + nPilot
    if (mod(i, pilotStep+1) == 0)
        cPilot(i) = 5;%constSymb(1);
    else
        cPilot(i) = c(j);
        j = j + 1;
    end
end

theta = zeros(length(c), 1);
theta(1) = rand()*2*pi;
for i = 2:nSymb
    theta(i) = mod(theta(i-1) + randn() * sqrt(phasevar),2*pi);
end

r = c .* exp(theta*1i) + randn(nSymb,1) * sqrt(noisevar) + randn(nSymb,1) * sqrt(noisevar) * 1i;

Pd = zeros(nSymb,M) + 1/M;

alpha = Pd * constSymb;
beta = Pd.^2 * constSymb;

a_f = zeros(nSymb,1);
for i = 2:nSymb
    a_f(i) = a_f(i-1) + (2 * r(i - 1) * conj(alpha(i - 1))) / (2 * noisevar + beta(i - 1) - abs(alpha(i - 1)))^2;
end

a_b = zeros(nSymb,1);
for i = nSymb-1:-1:1
    a_b(i) = a_b(i + 1) + 2 * (r(i + 1) * conj(alpha(i + 1))) / ((2 * noisevar + beta(i + 1) - abs(alpha(i + 1))^2));
end

Pu = zeros(nSymb, M);
for i = 1:nSymb
    for j = 1:M
        Pu(i,j) = exp(-abs(constSymb(j))^2 / (2 * noisevar)) * besseli(0, abs(a_f(i) + a_b(i) + (r(i) * conj(constSymb(j)))/noisevar));
    end
end

mu = zeros(2, length(encodedData)) + 0.5;

% Each possible combination of bits composing c_k
comb = zeros(2^(log2(M)-1), log2(M)-1);
for i = 0:2^(log2(M)-1)-1
    comb(i+1,:) = double(dec2bin(i, log2(M)-1)) - 48;
end
sizeComb = size(comb);

% Higher part of the FG: from Pu(c_k) to bit LLR 
bitCounter = 1;
newMu = zeros(1,length(encodedData));
for i = 1:nSymb
    for position = 1:log2(M) % Each possible position of bit b_k
        bitProb = zeros(1,2);
        for b_k = [0,1] % Bit b_k composing c_k
            for j = 1:sizeComb(1)
                % From bits to integer representation taking into accout the position of bit c_k
                exponent = sizeComb(2)+1;
                int = 0;
                curs = 1;
                for jj = 1:sizeComb(2)+1
                    if(jj == position)
                        int = int + b_k * 2^(exponent-jj);
                    else
                        int = int + comb(j,curs) * 2^(exponent-jj);
                        curs = curs + 1;
                    end
                end
                for m = 1:M % Each possible c_k
                    if (constSymb(m) == constSymbMod(int+1)) % Indicator function
                        for jj = 1:sizeComb(2)
                            bitProb(b_k+1) = bitProb(b_k+1) + Pu(i,m) * mu(b_k + 1, mod(i + position + jj, i*M)+1);
                        end
                    end
                end
            end
        end
        newMu(bitCounter) = bitProb(1) / bitProb(2);
        bitCounter = bitCounter + 1;
    end
end