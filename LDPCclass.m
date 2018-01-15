% Using handle is possible to modify the class properties within the class, and the changes are also seen in the main workspace
classdef LDPCclass < handle
    
    properties
        ParityCheckMatrix         % The Parity Check matrix of the LDPC code
        MaxNumberIterations       % The maximum number of iterations
        matCNI = 'matCNI.mat'     % Matrix containing check nodes indeces
        matVNI = 'matVNI.mat'     % Matrix containing variable nodes indeces
        MessageCheckToVariable    % Message from a check node to a variable node
        MessageVariableToCheck    % Message from a variable node to a check node
        nVN                       % Number of variable nodes
        nCN                       % Number of check nodes
        indexCN                   % Indeces of check nodes connected to each variable node (equivalent to CheckNodesIndex, but in matrix form, more efficient)
        indexVN                   % Indices of variable nodes connected to each check node (equivalent to VariableNodesIndex, but in matrix form, more efficient)
        Syndrome                  % LDPC Syndrome
        nCNmaxSize                % Maximum number of CN linked to a VN
        nVNmaxSize                % Maximum number of VN linked to a CN
        CheckNodesIndexLength     % Number of CN linked to a VN (for each VN)
        VariableNodesIndexLength  % Number of VN linked to a CN (for each CN)
    end
    
    methods
        %%%%%%% LDPC decoding function %%%%%%%
        % 1) Input: the object referred to the class; the bit LLR vector
        % 2) Output: the bit LLR calculated after the belief propagation of the graph; the weighted syndrome
        function [softOutput, wsyn]  = LDPC_decode(obj, bitLLR)
            % Save a local copy of variables of the object to improve speed performance
            ParityCheckMatrix_local = obj.ParityCheckMatrix;
            MessageVariableToCheck_local = obj.MessageVariableToCheck;
            MessageCheckToVariable_local = obj.MessageCheckToVariable;
            indexCN_local = obj.indexCN;
            indexVN_local = obj.indexVN;
            CheckNodesIndexLength_local = obj.CheckNodesIndexLength;
            VariableNodesIndexLength_local = obj.VariableNodesIndexLength;
            
            % NOTE: numerical problems related to the fact that MessageVariableToCheck martix may contains zeros, if there one variable
            % node has less connections with respect to the others. In this case, in order to make the algorithm usable in all cases,
            % we set zero elements in the matrix equal to a large number, e.g. 50, such that the tanh(50 / 2) = 1, and the one is irrelevant
            % in the product. Thanks to this, it is possible to extend the dimension of MessageVariableToCheck to embrace all possible
            % cases, wrt to different parity check matrices, without any problem.
            MessageVariableToCheck_local(MessageVariableToCheck_local == 0) = 50; % Make the tanh = 1 for the product
            
            % For loop on the number of iterations
            for numIt = 1:obj.MaxNumberIterations
                
                % Check Nodes Update
                count = ones(obj.nCN, 1);
                sumTot = sum(MessageCheckToVariable_local, 2); % Sum of incoming messages (from check nodes) for each variable node
                sumToCN = - MessageCheckToVariable_local + bitLLR + sumTot;
                % count1 = ones(1, obj.nCN);
                height = length(MessageVariableToCheck_local);
                for i = 1:obj.nVN % For each variable node
                    for j = 1:CheckNodesIndexLength_local(i) % For each check node connected to the variable node of index i
                        % The message from the variable node to the check node is equal to the sum of messages to the check node minus the
                        % message from the considered check node (of index j), plus the bitLLR from the detector (input of the decoder).
                        % Each row of the matrix MessageVariableToCheck represents a check node, and columns are the variable nodes
                        % connected to: the (i,j) element represent the message from variable node j to check node i.
                        MessageVariableToCheck_local(indexCN_local(i,j) + height * (count(indexCN_local(i,j)) - 1)) = sumToCN(i,j); % Select element by index to improve performance
                        % Update the indices: each check node receives a message from different variable nodes, if each column represents
                        % a variable node, once a message is received and stored in the (i,1) cell of the matrix, the message from the next
                        % variable node (to the same check node), has to be located in position (i,j+1).
                        count(indexCN_local(i,j)) = count(indexCN_local(i,j)) + 1;
                    end
                    % Try to vectorize the inner for loop, poor performance
                    % v = 1:CheckNodesIndexLength_local(i);
                    % % ind = sub2ind(size(MessageVariableToCheck_local), indexCN_local(i,1:nCNforVN), count1(indexCN_local(i,1:nCNforVN)));
                    % ind = indexCN_local(i,v) + height(1) * (count1(indexCN_local(i,v)) - 1);
                    % MessageVariableToCheck_local(ind) = sumToCN(i,v);
                    % count1(indexCN_local(i,v)) = count1(indexCN_local(i,v)) + 1;
                end
                
                % Variable Nodes Update
                count = ones(1,obj.nVN);
                prodTot = prod(tanh(MessageVariableToCheck_local / 2), 2);
                tanhMat = tanh(MessageVariableToCheck_local / 2);
                atanhMat = 2 * atanh((tanhMat ./ prodTot).^(-1));
                atanhMat(atanhMat == Inf) = 18.715;
                atanhMat(atanhMat == -Inf) = -18.715;
                height = length(MessageCheckToVariable_local);
                for i = 1:obj.nCN
                    for j = 1:VariableNodesIndexLength_local(i)
                        MessageCheckToVariable_local(indexVN_local(i,j) + height * (count(indexVN_local(i,j)) - 1)) = atanhMat(i,j); % Select element by index to improve performance
                        count(indexVN_local(i,j)) = count(indexVN_local(i,j)) + 1;
                    end
                end
                
                % Final APP soft information
                softOutput = bitLLR + sum(MessageCheckToVariable_local, 2);
                
                % Final APP hard information (bits)
                hardBits = softOutput;
                hardBits(hardBits > 0) = 0;
                hardBits(hardBits < 0) = 1;
                
                % Check Syndrome
                wsyn = sum(mod(ParityCheckMatrix_local * hardBits, 2));
                % If numerical problems arise, stop decoding. In main, a check on wsyn = Nan should be performed
                if(isnan(wsyn))
                    break;
                end
                % If the weighted syndrome is equel to 0, the correct codeword has been found, hence stop decoding
                if(wsyn == 0)
                    break;
                end
            end
            % At the end save the state of the decider
            obj.MessageVariableToCheck = MessageVariableToCheck_local;
            obj.MessageCheckToVariable = MessageCheckToVariable_local;
        end
        
        % Initialize the LDPC class
        % The fucntion takes as input the LDPC matrix (sparse logic) and produces
        % as output two files, specified after, in form of cells array.
        % matCn is the filename containing the indices of the check nodes connected to
        % a specific variable node (for each variable node).
        % matVn is the filename containing the indices of the variable nodes connected
        % to a specific check node (for each check node).
        % Take indices of nonzero elements of ldpc matrix,  save them into a cell array.
        function LDPC_initialization(obj, p, It)
            
            obj.ParityCheckMatrix = p; % Parity Check Matrix
            
            obj.MaxNumberIterations = It; % Maximum number of decoding iterations
            
            matDim = size(obj.ParityCheckMatrix); % Dimension of the Parity Check Matrix
            obj.nVN = matDim(2); % Number of variable nodes
            obj.nCN = matDim(1); % NUmber of check nodes
            
            % Save indices of check nodes for each variable node
            nVNmaxSize_local = 0;
            CN = cell(1, obj.nVN);
            CheckNodesIndexLength_local = zeros(obj.nVN, 1);
            for i = 1:obj.nVN
                CN{i} = find(obj.ParityCheckMatrix(:,i));
                if (length(CN{i}) > nVNmaxSize_local)
                    nVNmaxSize_local = length(CN{i});
                end
                CheckNodesIndexLength_local(i) = length(CN{i});
            end
            
            % Save indices of variable nodes for each check node
            ParityCheckMatrixT = transpose(obj.ParityCheckMatrix);
            nCNmaxSize_local = 0;
            VN = cell(1, obj.nCN);
            VariableNodesIndexLength_local = zeros(obj.nCN, 1);
            for i = 1:obj.nCN
                VN{i} = find(ParityCheckMatrixT(:,i));
                if (length(VN{i}) > nCNmaxSize_local)
                    nCNmaxSize_local = length(VN{i});
                end
                VariableNodesIndexLength_local(i) = length(VN{i});
            end
            
            obj.MessageCheckToVariable = zeros(obj.nVN, nVNmaxSize_local); % Message from check to variable nodes
            obj.MessageVariableToCheck = zeros(obj.nCN, nCNmaxSize_local); % Message from variable to check nodes
            
            % Indeces of check nodes saved in a matrix, to improve speed performance
            indexCN_local = zeros(obj.nVN, nVNmaxSize_local);
            for i = 1:obj.nVN
                indexCN_local(i,1:length(CN{i})) = CN{i};
            end
            % Indeces of variable nodes saved in a matrix, to improve speed performance
            indexVN_local = zeros(obj.nVN, nCNmaxSize_local);
            for i = 1:obj.nCN
                indexVN_local(i,1:length(VN{i})) = VN{i};
            end
            
            obj.CheckNodesIndexLength = CheckNodesIndexLength_local;
            obj.VariableNodesIndexLength = VariableNodesIndexLength_local;
            
            obj.nCNmaxSize = nCNmaxSize_local;
            obj.nVNmaxSize = nVNmaxSize_local;
            
            obj.indexCN = indexCN_local;
            obj.indexVN = indexVN_local;
        end
        
        % Reset internal messages of the LDPC graph. Must be called each time a new decoding stage starts. If joint detection/decoding
        % is performed, this method must be called before the iterations of the joint detector/decoder.
        function LDPC_reset(obj)
            obj.MessageCheckToVariable = zeros(obj.nVN, obj.nVNmaxSize); % Message from check to variable nodes
            obj.MessageVariableToCheck = zeros(obj.nCN, obj.nCNmaxSize); % Message from variable to check nodes
        end
    end
    
end