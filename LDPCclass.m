% Using handle is possible to modify the class properties within the class, and the changes are also seen in the main workspace
classdef LDPCclass < handle
    
    properties
        ParityCheckMatrix       % The Parity Check matrix of the LDPC code
        MaxNumberIterations     % The maximum number of iterations
        CheckNodesIndex         % Indeces of check nodes connected to each variable node
        matCNI = 'matCNI.mat'   % Matrix containing check nodes indeces
        VariableNodesIndex      % Indices of variable nodes connected to each check node
        matVNI = 'matVNI.mat'   % Matrix containing variable nodes indeces
        MessageCheckToVariable  % Message from a check node to a variable node
        MessageVariableToCheck  % Message from a variable node to a check node
        nVN                     % Number of variable nodes
        nCN                     % Number of check nodes
        indexCN                 % Indeces of check nodes connected to each variable node (equivalent to CheckNodesIndex, but in matrix form, more efficient)
        indexVN                 % Indices of variable nodes connected to each check node (equivalent to VariableNodesIndex, but in matrix form, more efficient)
        Syndrome                % LDPC Syndrome
    end
    
    methods
        % The fucntion takes as input the LDPC matrix (sparse logic) and produces
        % as output two files, specified after, in form of cells array.
        % matCn is the filename containing the indices of the check nodes connected to
        % a specific variable node (for each variable node).
        % matVn is the filename containing the indices of the variable nodes connected
        % to a specific check node (for each check node).
        % Take indices of nonzero elements of ldpc matrix,  save them into a cell
        % array and then the cell array into a file.
        % The string select if the function has to save or read the matrices.
        function LDPC_saveNodesConnections(obj)
            % Save indices of check nodes for each variable node
            CN = cell(1, obj.nVN);
            for i = 1:obj.nVN
                CN{i} = find(obj.ParityCheckMatrix(:,i))';
            end
            save(obj.matCNI, 'CN')
            
            % Save indices of variable nodes for each check node
            VN = cell(1, obj.nCN);
            for i = 1:obj.nCN
                VN{i} = find(obj.ParityCheckMatrix(i,:));
            end
            save(obj.matVNI, 'VN')
        end
        
        function LDPC_loadNodesConnections(obj)
            load(obj.matCNI);
            load(obj.matVNI);
            obj.CheckNodesIndex = CN;
            obj.VariableNodesIndex = VN;
        end
        
        function [softOutput, wsyn]  = LDPC_decode(obj, bitLLR)
            % Save a local copy of variables of the object to improve speed performance
            ParityCheckMatrix_local = obj.ParityCheckMatrix;
            MessageVariableToCheck_local = obj.MessageVariableToCheck;
            MessageCheckToVariable_local = obj.MessageCheckToVariable;
            CheckNodesIndex_local = obj.CheckNodesIndex;
            VariableNodesIndex_local = obj.VariableNodesIndex;
            indexCN_local = obj.indexCN;
            indexVN_local = obj.indexVN;
            
            % For loop on the number of iterations
            for numIt = 1:obj.MaxNumberIterations
                
                % Check Nodes Update
                count = ones(1, obj.nCN);
                % count1 = ones(1, obj.nCN);
                sumTot = sum(MessageCheckToVariable_local, 2); % Sum of incoming messages (from check nodes) for each variable node
                % height = size(MessageVariableToCheck_local);
                for i = 1:obj.nVN % For each variable node
                    for j = 1:length(CheckNodesIndex_local{i}) % For each check node connected to the variable node of index i
                        % The message from the variable node to the check node is equal to the sum of messages to the check node minus the
                        % message from the considered check node (of index j), plus the bitLLR from the detector (input of the decoder).
                        % Each row of the matrix MessageVariableToCheck represents a check node, and columns are the variable nodes
                        % connected to: the (i,j) element represent the message from variable node j to check node i.
                        MessageVariableToCheck_local(indexCN_local(i,j), count(indexCN_local(i,j))) = bitLLR(i) + sumTot(i) - MessageCheckToVariable_local(i,j);
                        % Update the indices: each check node receives a message from different variable nodes, if each column represents
                        % a variable node, once a message is received and stored in the (i,1) cell of the matrix, the message from the next
                        % variable node (to the same check node), has to be located in position (i,j+1).
                        count(indexCN_local(i,j)) = count(indexCN_local(i,j)) + 1;
                    end
                    % Try to vectorize the inner for
%                     nCNforVN = length(CheckNodesIndex_local{i});
%                     % ind = sub2ind(size(MessageVariableToCheck_local), indexCN_local(i,1:nCNforVN), count1(indexCN_local(i,1:nCNforVN)));
%                     ind = indexCN_local(i,1:nCNforVN) + height(1) .* (count1(indexCN_local(i,1:nCNforVN)) - 1);
%                     MessageVariableToCheck_local(ind) = bitLLR(i) + sumTot(i) - MessageCheckToVariable_local(i,1:nCNforVN);
%                     count1(indexCN_local(i,1:nCNforVN)) = count1(indexCN_local(i,1:nCNforVN)) + 1;
                end
                
                % Variable Nodes Update
                count = ones(1,obj.nVN);
                prodTot = prod(tanh(MessageVariableToCheck_local / 2), 2);
                for i = 1:obj.nCN
                    for j = 1:length(VariableNodesIndex_local{i})
                        MessageCheckToVariable_local(indexVN_local(i,j), count(indexVN_local(i,j))) = 2 * atanh(prodTot(i) / tanh(MessageVariableToCheck_local(i,j) / 2));
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
                check = mod(ParityCheckMatrix_local * hardBits, 2);
                wsyn = sum(check);
                % Is numerical problems arise, stop decoding. In main, a check on wsyn = Nan should be performed
                if(isnan(wsyn))
                    break;
                end
                % If the weighted syndrome is equel to 0, the correct codeword has been found, hence stop decoding
                if(sum(check) == 0)
                    break;
                end
            end
            % At the end save the state of the decider
            obj.MessageVariableToCheck = MessageVariableToCheck_local;
            obj.MessageCheckToVariable = MessageCheckToVariable_local;
        end
        
        % Initialize the LDPC class
        function LDPC_initialization(obj, p, It, cmd)
            
            obj.ParityCheckMatrix = p; % Parity Check Matrix
            
            obj.MaxNumberIterations = It; % Maximum number of decoding iterations
            
            matDim = size(obj.ParityCheckMatrix); % Dimension of the Parity Check Matrix
            obj.nVN = matDim(2); % Number of variable nodes
            obj.nCN = matDim(1); % NUmber of check nodes
            
            if (cmd == 'save')
                % Save indices of check nodes for each variable node
                CN = cell(1, obj.nVN);
                for i = 1:obj.nVN
                    CN{i} = find(obj.ParityCheckMatrix(:,i))';
                end
                save(obj.matCNI, 'CN')
                
                % Save indices of variable nodes for each check node
                VN = cell(1, obj.nCN);
                for i = 1:obj.nCN
                    VN{i} = find(obj.ParityCheckMatrix(i,:));
                end
                save(obj.matVNI, 'VN')
                
                obj.CheckNodesIndex = CN;
                obj.VariableNodesIndex = VN;
            end
            
            if (cmd == 'load')
                load(obj.matCNI);
                load(obj.matVNI);
                obj.CheckNodesIndex = CN;
                obj.VariableNodesIndex = VN;
            end
            
            obj.MessageCheckToVariable = zeros(obj.nVN, 8); % Message from check to variable nodes
            obj.MessageVariableToCheck = zeros(obj.nCN, 7); % Message from variable to check nodes
            
            CheckNodesIndex_local = obj.CheckNodesIndex; % Indeces of check nodes saved in a cell array
            VariableNodesIndex_local = obj.VariableNodesIndex; % Indeces of variable nodes saved in a cell array
            
            % Indeces of check nodes saved in a matrix, to improve speed performance
            indexCN_local = zeros(obj.nVN, 30);
            for i = 1:obj.nVN
                indexCN_local(i,1:length(CheckNodesIndex_local{i})) = CheckNodesIndex_local{i};
            end
            % Indeces of variable nodes saved in a matrix, to improve speed performance
            indexVN_local = zeros(obj.nVN, 30);
            for i = 1:obj.nCN
                indexVN_local(i,1:length(VariableNodesIndex_local{i})) = VariableNodesIndex_local{i};
            end
            
            obj.indexCN = indexCN_local;
            obj.indexVN = indexVN_local;
        end
        
        function LDPC_reset(obj)
            obj.MessageCheckToVariable = zeros(obj.nVN, 8); % Message from check to variable nodes
            obj.MessageVariableToCheck = zeros(obj.nCN, 7); % Message from variable to check nodes
        end
    end
    
end
