clc
clear

rng(15);

%%%%%%%%%%%%%%% System Parameters %%%%%%%%%%%%%%%
M = 4;                                                          % Constellation Cardinality
R = 1/2;                                                        % Code Rate
p = dvbs2ldpc(R);                                               % LDPC parity check matrix
hEnc = comm.LDPCEncoder(p);                                     % Definition of the LDPC encoder
hMod = comm.PSKModulator(M, 'BitInput',true, 'PhaseOffset', 0); % Definition of the Modulator
obj = LDPCclass;                                                % Definition of the LDPC class
% LDPC_initialization(obj, p, 1, 'save');                       % Initialize LDPC class ('save' to save matrices, must be called once)
LDPC_initialization(obj, p, 1, 'load');                         % Initialize LDPC class ('load' to load matrices saved before)
codewordLength = 64800;                                         % Number of coded bits
messageLength = 64800 * R;                                      % Length of the data message
nSymb = codewordLength / log2(M);                               % Number of Symbols
maxIt = 100;

%%%%%%%%%%% Constellation Management %%%%%%%%%%%%
constSymb = constellation(hMod);                                % Constellation symbols
constBits = zeros(M, log2(M));                                  % Constellation bit mapping
constSymbMod = zeros(M, log2(M));                               % Modified constellation bit mapping (wrt to integer representation)
for i = 0:2^(log2(M))-1
    symbol = hMod(double(dec2bin(i, log2(M))).' - 48);
    index = find(constSymb == symbol);
    constBits(index,:) = double(dec2bin(i, log2(M))) - 48;
    constSymbMod(i+1) = symbol;
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
EbN0 = 6:6;
for l = 1:length(EbN0)
    
    noisevar = 1 / (log2(M) * R * 10^(EbN0(l)/10)); % Noise variance
    
    trials = 0;
    errors = 0;
    while(errors < 100)
        
        %%%%%%%%%% Information Part %%%%%%%%%%
        data = logical(randi([0 1], messageLength, 1)); % Message bits
        encodedData = hEnc(data);                       % Encoded bits
        c = hMod(encodedData);                          % Modulated symbols
        
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
        
        LDPC_reset(obj); % Reset the state of the LDPC
        
        %%%%%%%%% Tikhonov Algorithm %%%%%%%%%
        for numIt = 1:maxIt
            
            bitProbT = transpose(bitProb);     % Transpose bit probability matrix since comm.toolbox objects need column vectors
            
            alpha = Pd * constSymb;            % First order moment of ck distribution
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
            bitCounter = 1;
            bitLLR = zeros(1,length(encodedData));
            exponent = sizeComb(2) + 1;
            probBit = zeros(1,2);
            for i = 0:length(encodedData) - 1
                % positions = log2(M) * floor(i/log2(M)) + 1:log2(M) * floor(i/log2(M)) + 1 + log2(M) - 1;
                % positions(mod(i, log2(M)) + 1) = [];
                probBit = probBit - probBit; % Equivalent to probBit = zeros(1,2) but faster
                for bk = 0:1 % Bit bk composing ck, want to calculate the LLR
                    for j = 1:sizeComb(1)
                        % From bits to integer representation taking into accout the position of bit ck
                        int = 0;
                        curs = 1;
                        for jj = 1:sizeComb(2) + 1
                            if(jj == mod(i, log2(M)) + 1)
                                int = int + bk * 2^(exponent-jj);
                            else
                                int = int + comb(j, curs) * 2^(exponent-jj);
                                curs = curs + 1;
                            end
                        end
                        for m = 1:M % Each possible c_k
                            if (constSymb(m) == constSymbMod(int + 1)) % Indicator function
                                mu = 1;
                                for jj = 1:sizeComb(2)
                                    mu = mu * bitProbT(comb(j, jj) + 1, positionsMat(i + 1, jj));
                                end
                                probBit(bk + 1) = probBit(bk + 1) + mu * pilotlessPu(floor(i/log2(M)) + 1, m);
                            end
                        end
                    end
                end
                bitLLR(bitCounter) = log(probBit(1) / probBit(2)); % Calculate the bit LLR
                bitCounter = bitCounter + 1;
            end
            
            [receivedLLR, wsyn] = LDPC_decode(obj, bitLLR');
            
            %%% From bit LLR to bit probabilities %%%
            bitProb(:,1) = exp(receivedLLR) ./ (1 + exp(receivedLLR));  % Probability of bit to be = 0
            bitProb(:,2) = 1 ./ (1 + exp(receivedLLR));                 % Probability of bit to be = 1
            
            %%% From bit to symbol probabilities %%%
            pilotlessPd = ones(nSymbPilot, M);
            ii = 1;
            for i = 1:log2(M):codewordLength
                if (mod(ii, pilotStep + 1) == 0)
                    ii = ii + 1;
                end
                for j = 1:M
                    for z = 1:log2(M)
                        pilotlessPd(ii, j) = pilotlessPd(ii, j) * bitProb(i + z - 1, constBits(j,z) + 1);
                    end
                end
                ii = ii + 1;
            end
            
            %%% Pilot Insertion %%%
            Pd(pilotStep + 1:pilotStep + 1:end, :) = 0;           % Not pilot symbols (at pilot positions) have probability 0
            Pd(pilotStep + 1:pilotStep + 1:end, pilotIndex) = 1;  % Pilot symbols (at pilot positions) have probability 1
            
            if(wsyn == 0 || isnan(wsyn))
                break;
            end
            fprintf('wsyn: %d \n', wsyn);
        end
        
        %%%%%%%%%% Error Counting %%%%%%%%%%
        if(~isnan(wsyn))
            receivedBits = receivedLLR(1:length(data));
            receivedBits(receivedBits > 0) = 0;
            receivedBits(receivedBits < 0) = 1;
            errors = sum(mod(receivedBits + data, 2));
            trials = trials + messageLength;
        end
        fprintf('EbN0: %f \t BER: %f \n', EbN0(l), errors/trials);
    end
    fprintf('EbN0: %f \t BER: %f \n', EbN0(l), errors/trials);
end
