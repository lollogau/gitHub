% The function takes as input:
% - the bit LLR vector from which it calculates the symbol probabilities
% - the constellation bit mapping: a matrix in which each row contains the
%   bits (one per each cell) of each symbol of the used constellation.
%   The mapping rule, i.e., Gray or another one, is implicitly specified in
%   the matrix.
% The output is a matrix. Each row i of the matrix contains the symbol probabilities
% related to the i-th bit LLR. The column j specifies the j-th
% constellation symbol considered. 

function symbProb = bitLLR2symbProb(bitLLR, constellationBits)

    %%% Constellation Order %%%
    M = length(constellationBits);
    
    %%% From bit LLR to bit probabilities %%%
    bitProb(:,1) = exp(bitLLR) ./ (1 + exp(bitLLR));     % Probability of bit to be = 0
    bitProb(:,2) = 1 ./ (1 + exp(bitLLR));               % Probability of bit to be = 1

    %%% From bit LLR to symbol probabilities %%%
    symbProb = ones(length(bitLLR)/log2(M), length(constellationBits));
    ii = 1;
    for i = 1:log2(M):length(bitLLR)
        for j = 1:M
            for z = 1:log2(M)
                symbProb(ii, j) = symbProb(ii, j) * bitProb(i + z - 1, constellationBits(j,z) + 1); % Implicit indicator function given by the selection of the j-th row of constBits
            end
        end
        ii = ii + 1;
    end
    
end

