clc
clear

rng(15);                                                        % Noise generator
fid = fopen('results_alpha_1.txt','w');

%%%%%%%%%%%%%%% System Parameters %%%%%%%%%%%%%%%
M = 4;                                                          % Constellation Cardinality
R = 1/2;                                                        % Code Rate
p = dvbs2ldpc(R);                                               % LDPC parity check matrix
hEnc = comm.LDPCEncoder(p);                                     % Definition of the LDPC encoder
hMod = comm.PSKModulator(M, 'BitInput',true, 'PhaseOffset', 0); % Definition of the Modulator
obj = LDPCclass;                                                % Definition of the LDPC class
LDPC_initialization(obj, p, 1);                                 % Initialize LDPC class
codewordLength = 64800;                                         % Number of coded bits
messageLength = 64800 * R;                                      % Length of the data message
nSymb = codewordLength / log2(M);                               % Number of Symbols
maxIt = 100;                                                    % Number of detector/decoder iterations

%%%%%%%%%%% Constellation Management %%%%%%%%%%%%
constSymb = constellation(hMod);                                % Constellation symbols
constBits = zeros(M, log2(M));                                  % Constellation bit mapping
constIntToGray = zeros(M, 1);
for i = 0:2^(log2(M))-1
    symbol = hMod(double(dec2bin(i, log2(M))).' - 48);
    index = find(constSymb == symbol);
    constBits(index,:) = double(dec2bin(i, log2(M))) - 48;
    constIntToGray(i + 1) = index;                              % Modified constellation bit mapping (wrt to integer representation). Indicates the position of the integer number, represented by the index of the vector
end
release(hMod);

comb = zeros(2^(log2(M)-1), log2(M)-1);                         % Each possible combination of bits composing ck
for i = 0:2^(log2(M)-1)-1
    comb(i+1,:) = double(dec2bin(i, log2(M)-1)) - 48;
end
sizeComb = size(comb);

positionsMat = zeros(codewordLength, log2(M) - 1);
for i = 0:codewordLength - 1
    positions = log2(M) * floor(i/log2(M)) + 1:log2(M) * floor(i/log2(M)) + 1 + log2(M) - 1;
    positions(mod(i, log2(M)) + 1) = [];
    positionsMat(i + 1, :) = positions;
end

exponents = zeros(sizeComb(2), log2(M));
for i = 0:log2(M) - 1
    bitPosition = mod(i, log2(M)) + 1; % Position of bit bk
    exponents(:, i + 1) = [1:bitPosition-1, bitPosition+1:log2(M)];
end

%%%%%%%%%%%%%%%% Pilot Symbols %%%%%%%%%%%%%%%%%
pilotIndex = 1;                                                 % The pilot symbol is the symbol of the constellation in position pilotIndex of constSymb
pilotSymbol = constSymb(pilotIndex);                            % Pilot symbol
pilotStep = 10;                                                 % Step between pilot symbols
nPilot = nSymb / pilotStep;                                     % Number of pilot symbols
nSymbPilot = nSymb + nPilot;                                    % NUmber of symbols + pilot symbols

%%%%%%%%%%%% Phase Noise Parameters %%%%%%%%%%%%
sigmaDelta = 0.3;
phasevar = (sigmaDelta / 180 * pi)^2;                           % Variance of the phase in radiant

%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%
EbN0 = 3:0.1:4;
for l = 1:length(EbN0)
    
    noisevar = 1 / (log2(M) * R * 10^(EbN0(l)/10));             % Noise variance
    
    trials = 0;
    errors = 0;
    pckErrors = 0;
    while(pckErrors < 100)
        
        %%%%%%%%%% Information Part %%%%%%%%%%
        data = logical(randi([0 1], messageLength, 1));         % Message bits
        encodedData = hEnc(data);                               % Encoded bits
        c = hMod(encodedData);                                  % Modulated symbols
        
        %%%%%%% Pilot symbols insertion %%%%%%
        cPilot = c;
        cPilot = reshape(cPilot, pilotStep, nSymb/pilotStep);
        cPilot(pilotStep + 1,:) = pilotSymbol;
        cPilot = reshape(cPilot, nSymbPilot, 1);
        
        %%%%%%%%% Phase noise process %%%%%%%%
        theta = zeros(length(cPilot), 1);
        theta(1) = rand()*2*pi;
        for i = 2:nSymbPilot
            theta(i) = mod(theta(i-1) + randn() * sqrt(phasevar),2*pi);
        end
        
        %%%%%%%%%%%%%% Channel %%%%%%%%%%%%%%%
        r = cPilot .* exp(theta*1i) + randn(nSymbPilot,1) * sqrt(noisevar) + randn(nSymbPilot,1) * sqrt(noisevar) * 1i;
        % r = cPilot + randn(nSymbPilot,1) * sqrt(noisevar) + randn(nSymbPilot,1) * sqrt(noisevar) * 1i;  

        %%%%%%%%% Symbol probability %%%%%%%%%
        Pd = zeros(nSymbPilot, M) + 1/M;                      % All symbols are equally probable at first iteration
        Pd(pilotStep + 1:pilotStep + 1:end, :) = 0;           % Not pilot symbols (at pilot positions) have probability 0
        Pd(pilotStep + 1:pilotStep + 1:end, pilotIndex) = 1;  % Pilot symbols (at pilot positions) have probability 1
        
        %%%%%%%%%% Bits probability %%%%%%%%%%
        bitProb = zeros(length(encodedData), 2) + 0.5;        % At first iteration all bits are equally probable
        
        LDPC_reset(obj);                                      % Reset the state of the LDPC
        
        bitLLR = zeros(codewordLength, 1);                    % Bit LLR vector
        exponent = sizeComb(2) + 1;                           % Maximum exponent when converting from bit to integer notation + 1
        probBit = zeros(1,2);                                 % Bit probability vector, to store and then calculate the LLR
        
        %%%%%%%%% Tikhonov Algorithm %%%%%%%%%
        for numIt = 1:maxIt
            
            bitProbT = transpose(bitProb);     % Transpose bit probability matrix since comm.toolbox objects need column vectors
            
            alpha = Pd * constSymb;            % First order moment of ck distribution
            % [maxPd, I] = max(Pd,[],2);
            % alpha = constSymb(I);
            beta = Pd * abs(constSymb).^2;     % Second order moment of ck distribution
            
            %%% Forward recursion (32) %%%
            aForw = zeros(nSymbPilot,1);
            for i = 2:nSymbPilot
                aForw(i) = aForw(i-1) + (2 * r(i - 1) * conj(alpha(i - 1))) / (2 * noisevar + beta(i - 1) - abs(alpha(i - 1))^2);
                aForw(i) = aForw(i) / (1 + phasevar * abs(aForw(i)));
            end
            
            %%% Backward recursion (34) %%%
            aBack = zeros(nSymbPilot,1);
            for i = nSymbPilot-1:-1:1
                aBack(i) = aBack(i + 1) + 2 * (r(i + 1) * conj(alpha(i + 1))) / (2 * noisevar + beta(i + 1) - abs(alpha(i + 1))^2);
                aBack(i) = aBack(i) / (1 + phasevar * abs(aBack(i)));
            end
            
            %%% Message to the LDPC code (35). In log domain approximate the Besseli with the argument %%%
            Pu = - (abs(transpose(constSymb)).^2 ./ (2 * noisevar)) + abs(aForw + aBack + (r * transpose(conj(constSymb)))/noisevar);
            Pu = exp(Pu - max(Pu,[],2)); % Normalization
            Pu = Pu ./ sum(Pu, 2);
            pilotlessPu = Pu;
            pilotlessPu(1+pilotStep:1+pilotStep:end,:) = [];
            
            %%% From soft symbols to bit LLR %%%
            for i = 0:codewordLength - 1
                bitPosition = mod(i, log2(M)) + 1;                         % Position of bit bk in the log2(M) bits
                probBit = probBit - probBit;                               % Equivalent to probBit = zeros(1,2) but faster
                for bk = 0:1                                               % Bit bk composing ck, want to calculate the LLR
                    for j = 1:sizeComb(1)
                        int = bk * 2^(exponent - bitPosition);             % Integer representation of log2(M) bits (first part)
                        for jj = 1:sizeComb(2)
                            int = int + comb(j, jj) * 2^(exponent - exponents(jj, bitPosition));   % Integer representation of log2(M) bits (second part)
                        end
                        
                        % Vectorized Version --------------------------
                        % int = bk * 2^(exponent-bitPosition) + comb(j, :)...
                        % * 2.^(exponent-exponents(:, bitPosition));
                        % ---------------------------------------------
                        
                        mu = 1;                                            % Message from the other nodes
                        for jj = 1:sizeComb(2)
                            mu = mu * bitProbT(comb(j, jj) + 1, positionsMat(i + 1, jj));
                        end
                        probBit(bk + 1) = probBit(bk + 1) + mu * pilotlessPu(floor(i/log2(M)) + 1, constIntToGray(int + 1)); % Bit probability
                    end
                end
                bitLLR(i + 1) = log(probBit(1) / probBit(2));              % Bit LLR
            end
            
            [receivedLLR, wsyn] = LDPC_decode(obj, bitLLR);
            
            %%% From bit LLR to bit probabilities %%%
            bitProb(:,1) = exp(receivedLLR) ./ (1 + exp(receivedLLR));     % Probability of bit to be = 0
            bitProb(:,2) = 1 ./ (1 + exp(receivedLLR));                    % Probability of bit to be = 1
            
            %%% From bit to symbol probabilities %%%
            Pd = ones(nSymbPilot, M);
            ii = 1;
            for i = 1:log2(M):codewordLength
                if (mod(ii, pilotStep + 1) == 0) % For pilot insertion
                    ii = ii + 1;
                end
                for j = 1:M
                    for z = 1:log2(M)
                        Pd(ii, j) = Pd(ii, j) * bitProb(i + z - 1, constBits(j,z) + 1); % Implicit indicator function given by the selection of the j-th row of constBits
                    end
                end
                ii = ii + 1;
            end
            
            %%% Pilot Insertion %%%
            Pd(pilotStep + 1:pilotStep + 1:end, :) = 0;           % Not pilot symbols (at pilot positions) have probability 0
            Pd(pilotStep + 1:pilotStep + 1:end, pilotIndex) = 1;  % Pilot symbols (at pilot positions) have probability 1
            
            % fprintf('wsyn: %d \n', wsyn);
            if(wsyn == 0 || isnan(wsyn))
                break;
            end
        end
        
        %%%%%%%%%% Error Counting %%%%%%%%%%
        if(~isnan(wsyn))
            receivedBits = receivedLLR(1:length(data));
            receivedBits(receivedBits > 0) = 0;
            receivedBits(receivedBits < 0) = 1;
            errors = errors + sum(mod(receivedBits + data, 2));
            trials = trials + messageLength;
            pckErrors = pckErrors + 1;
        end
        fprintf('EbN0: %f \t BER: %f \n', EbN0(l), errors/trials);
    end
    fprintf('EbN0: %f \t BER: %f \n', EbN0(l), errors/trials);
end
