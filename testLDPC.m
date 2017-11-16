clc
clear all

M = 4; % COnstellation Cardinality

hEnc = comm.LDPCEncoder;
hMod = comm.PSKModulator(M, 'BitInput',true, 'PhaseOffset', 0);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',5);
hDemod = comm.PSKDemodulator(4, 'BitOutput',true,'DecisionMethod','Approximate log-likelihood ratio','Variance', 1/10^(hChan.SNR/10));
hDec = comm.LDPCDecoder('OutputValue','Whole codeword','DecisionMethod','Soft decision');
hError = comm.ErrorRate;

const = constellation(hMod); % Constellation
constBits = zeros(M, log2(M)); % Vector containing the bits of constellation symbols
for i = 0:2^(log2(M))-1
    symbol = hMod(double(dec2bin(i, log2(M)))' - 48);
    index = find(const == symbol);
    constBits(index,:) = double(dec2bin(i, log2(M))) - 48;
    constSymbMod(i+1) = symbol;
end
release(hMod);

noisevar = 1 / 10^(hChan.SNR/10);
phasevar = 1 / 10^(5/10);

data = logical(randi([0 1], 32400, 1)); % Message bits
encodedData = hEnc(data); % Encoded bits
c = hMod(encodedData); % Modulated symbols
nSymb = length(c); % Number of Symbols

theta = zeros(length(c), 1);
theta(1) = rand()*2*pi;
for i = 2:nSymb
    theta(i) = mod(theta(i-1) + randn() * sqrt(phasevar),2*pi);
end

r = c .* exp(theta*1i) + randn(nSymb,1) * sqrt(noisevar);

Pd = zeros(nSymb,M) + 1/M;

alpha = Pd * const;
beta = Pd.^2 * const;

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
        Pu(i,j) = exp(-abs(const(j))^2 / (2 * noisevar)) * besseli(0, abs(a_f(i) + a_b(i) + (r(i) * const(j))/noisevar));
    end
end

mu = zeros(2, length(encodedData)) + 0.5;

% Each possible combination of bits composing c_k
comb = zeros(2^(log2(M)-1), log2(M)-1);
for i = 0:2^(log2(M)-1)-1
    comb(i+1,:) = double(dec2bin(i, log2(M)-1)) - 48;
end
sizeComb = size(comb);

c_k = 0;
position = 1;
for i = 1:nSymb
    for m = 1:M % Each possible c_k
        for j = 1:sizeComb(1)
            exponent = sizeComb(2)+1;
            int = 0;
            curs = 1;
            for jj = 1:sizeComb(2)+1
                if(jj == position)
                    int = int + c_k * 2^(exponent-jj);
                else
                    int = int + comb(j,curs) * 2^(exponent-jj);
                    curs = curs + 1;
                end
            end
            %v = insert(comb(j,:), c_k, position);
            %[~,index] = ismember(constBits, v, 'rows');
            %index = find (index == 1);
            %v = bin2dec(num2str(v));
            %int = 0;
            %for jj = 1:length(v)
                %int = int + v(jj) * 2^(length(v)-jj);
            %end
            if (constSymbMod(m) == int) % Indicator function
                newMu = 0;
                for jj = 1:length(comb(j))
                    newMu = newMu + Pu(i,m) * mu(mod(i + position +jj, i*M));
                end
                
            end
        end
    end
end

% function v = insert(v, n, nIndex)
%     v = [v(1:nIndex-1), n, v(nIndex:end)];
% end