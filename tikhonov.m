clc
clear all

M = 4; % COnstellation Cardinality
pilotIndex = 1; % The pilot symbol is the symbol of the constellation in position pilotIndex

hEnc = comm.LDPCEncoder;
hMod = comm.PSKModulator(M, 'BitInput',true, 'PhaseOffset', 0);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',2);
hDec = comm.LDPCDecoder('OutputValue','Whole codeword','DecisionMethod','Soft decision','MaximumIterationCount', 10);
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
phasevar = 1;

data = logical(randi([0 1], 32400, 1)); % Message bits
encodedData = hEnc(data); % Encoded bits
c = hMod(encodedData); % Modulated symbols
nSymb = length(c); % Number of Symbols

% Pilot symbols insertion
pilotSymbol = constSymb(pilotIndex);
pilotStep = 10;
nPilot = length(c) / pilotStep ; % Number of pilot symbols
cPilot = zeros(nSymb + nPilot, 1);
j = 1;
for i = 1:nSymb + nPilot
    if (mod(i, pilotStep+1) == 0)
        cPilot(i) = pilotSymbol;
    else
        cPilot(i) = c(j);
        j = j + 1;
    end
end
nSymb = length(cPilot); % The number of symbols takes into account the number of pilot symbols

theta = zeros(length(cPilot), 1);
theta(1) = rand()*2*pi;
for i = 2:nSymb
    theta(i) = mod(theta(i-1) + randn() * sqrt(phasevar),2*pi);
end

r = cPilot .* exp(theta*1i) + randn(nSymb,1) * sqrt(noisevar) + randn(nSymb,1) * sqrt(noisevar) * 1i;
% r = cPilot + randn(nSymb,1) * sqrt(noisevar) + randn(nSymb,1) * sqrt(noisevar) * 1i;

% Symbol probability
Pd = zeros(nSymb, M) + 1/M; % All symbols are equally probable at first iteration
Pd(pilotStep + 1:pilotStep + 1:end, :) = 0; % Other not pilot symbols (at pilot positions) have probability 0
Pd(pilotStep + 1:pilotStep + 1:end, pilotIndex) = 1; % The pilot symbols (at pilot positions) have probability 1

% Each possible combination of bits composing ck
comb = zeros(2^(log2(M)-1), log2(M)-1);
for i = 0:2^(log2(M)-1)-1
    comb(i+1,:) = double(dec2bin(i, log2(M)-1)) - 48;
end
sizeComb = size(comb);

bitProb = zeros(length(encodedData), 2) + 0.5; % At first iteration all bits are equally probable

for numIt = 1:50
    
    bitProbT = transpose(bitProb); % Transpose bit probability matrix since comm.toolbox objects need column vectors
    
    alpha = Pd * constSymb; % First order moment of ck distribution
    beta = Pd * abs(constSymb).^2; % Second order moment of ck distribution
    
    % Forward recursion
    aForw = zeros(nSymb,1);
    for i = 2:nSymb
        aForw(i) = aForw(i-1) + (2 * r(i - 1) * conj(alpha(i - 1))) / (2 * noisevar + beta(i - 1) - abs(alpha(i - 1))^2);
    end
    
    % Backward recursion
    aBack = zeros(nSymb,1);
    for i = nSymb-1:-1:1
        aBack(i) = aBack(i + 1) + 2 * (r(i + 1) * conj(alpha(i + 1))) / (2 * noisevar + beta(i + 1) - abs(alpha(i + 1))^2);
    end
    
    % Message to the LDPC code
    Pu = zeros(nSymb, M);
    for i = 1:nSymb
        for j = 1:M
            % Pu(i,j) = exp(-abs(constSymb(j))^2 / (2 * noisevar)) * besseli(0, abs(aForw(i) + aBack(i) + (r(i) * conj(constSymb(j)))/noisevar));
            Pu(i,j) = - (abs(constSymb(j))^2 / (2 * noisevar)) + abs(aForw(i) + aBack(i) + ((r(i) * conj(constSymb(j)))/noisevar));
        end
        % Pu(i,:) = Pu(i,:) / sum(Pu(i,:)); % Normalitazion
        Pu(i,:) = Pu(i,:) - max(Pu(i,:));
        Pu(i,:) = exp(Pu(i,:));
        Pu(i,:) = Pu(i,:) / sum(Pu(i,:));
    end
    
    % Higher part of the FG: from Pu(c_k) to bit LLR
    pilotlessPu = Pu;
    pilotlessPu(1+pilotStep:1+pilotStep:end,:) = [];
    bitCounter = 1;
    bitLLR = zeros(1,length(encodedData));
    for i = 1:length(encodedData)/log2(M)
        for position = 1:log2(M) % Each possible position of bit b_k
            probBit = zeros(1,2);
            for bk = [0,1] % Bit b_k composing c_k
                for j = 1:sizeComb(1)
                    % From bits to integer representation taking into accout the position of bit c_k
                    exponent = sizeComb(2)+1;
                    int = 0;
                    curs = 1;
                    for jj = 1:sizeComb(2)+1
                        if(jj == position)
                            int = int + bk * 2^(exponent-jj);
                        else
                            int = int + comb(j,curs) * 2^(exponent-jj);
                            curs = curs + 1;
                        end
                    end
                    for m = 1:M % Each possible c_k
                        if (constSymb(m) == constSymbMod(int+1)) % Indicator function
                            for jj = 1:sizeComb(2)
                                probBit(bk+1) = probBit(bk+1) + pilotlessPu(i,m) * bitProbT(bk + 1, mod(i + position + jj, i*M)+1);
                            end
                        end
                    end
                end
            end
            bitLLR(bitCounter) = log(probBit(1) / probBit(2)); % Calculate the bit LLR
            bitCounter = bitCounter + 1;
        end
    end
    
    receivedLLR = hDec(bitLLR');
    
    % From LLR to symbols probabilities
    bitProb(:,1) = exp(receivedLLR) ./ (1 + exp(receivedLLR)); % Probability of bit to be = 0
    bitProb(:,2) = 1 ./ (1 + exp(receivedLLR)); % Probability of bit to be = 1
    
    % Symbol probabilities
    Pd = ones(nSymb, M);
    ii = 1;
    for i = 1:log2(M):length(receivedLLR)
        if (mod(ii, pilotStep + 1) == 0)
            Pd(ii, :) = zeros();
            Pd(ii, pilotIndex) = pilotSymbol;
        else
            for j = 1:M
                for z = 1:log2(M)
                    Pd(ii, j) = Pd(ii, j) * bitProb(i + z - 1, constBits(j,z) + 1);
                end
            end
        end
        ii = ii + 1;
    end
    
    % Pd(pilotStep + 1:pilotStep + 1:end, :) = 0; % Other not pilot symbols (at pilot positions) have probability 0
    % Pd(pilotStep + 1:pilotStep + 1:end, pilotIndex) = 1; % The pilot symbols (at pilot positions) have probability 1
end
