function SINRPerRE = ComputePrecodedSINRPerRE(H,codebook,codebookIdxSetSizes,nVar,k,numLayers)
%

% Get the number of sc in BWP
csirsIndSubsLen = length(k);

% Allocate SINRPerRE
SINRPerRE = zeros([csirsIndSubsLen,numLayers,codebookIdxSetSizes]);

% Total number of element in codebook index set size
prodSetSizes = prod(codebookIdxSetSizes);

if csirsIndSubsLen > prodSetSizes
    % If the number of csirs subcarriers in BWP is greater than codebook 
    % index set size, consider the channel matrix as the page matrix
    for idx = 1:prodSetSizes
        % Check if all elements are zeros
        if any(codebook(:,:,idx),'all')
            % Compute linear SINR of all RE that contains CSI-RS per 
            % precoding matrix
            WTmp = codebook(:,:,idx);
            sinr = computePrecodedSINRSubroutine(H,WTmp,nVar);
            SINRPerRE(:,:,idx) = sinr;
        end
    end
else
    % If the number of csirs subcarriers in BWP is less than codebook 
    % index set size, consider the precoding matrix as the page matrix
    for idx = 1:csirsIndSubsLen
        % Compute the SINR (linear) values per CSI-RS RE for all
        % unrestricted precoding matrices
        HTmp = H(idx,:,:);
        sinr = computePrecodedSINRSubroutine(HTmp,codebook,nVar);
        % Assign sinr [k x l x m] (size of codebookIdxSetSizes into 
        % SINRPerRE slices because SINRPerRE(i,j,:,:,:) matches size(sinr) 
        % for each i,j
        SINRPerRE(idx,:) = sinr(:);
    end
end
end

%% Local Helper Fcn
function sinr = computePrecodedSINRSubroutine(H,W,nVar)
% Reference
% nrPrecodedSINR(H,nVar,W) from Mathworks
% Li, Ping, et al. "On the Distribution of SINR for the MMSE MIMO Receiver 
% and Performance Analysis." IEEE Transactions on Information Theory, 
% vol. 52, no. 1, 1 Jan. 2006, pp. 271–286, 
% https://doi.org/10.1109/tit.2005.860466. Accessed 23 Sept. 2023.

% Rerrange H matrix, so H is now nRxAnt-Pcsirs-K*L
H = permute(H,[2, 3, 1]); % Toss the BWP sc to the third dim to match W

R = pagemtimes(H,W);
[~,S,V] = pagesvd(R,'econ','vector');

SSqr = S.*S;
absVSqr = abs(V).^2;

% If H is 2D, compute SINR values using W as page
if size(H,3) == 1
    diagTerm = nVar./(nVar + SSqr + eps); % eps to prevent divide by 0
    diagTermPageT = pagetranspose(diagTerm);
    sumTerm = sum(absVSqr.*diagTermPageT,2);
    msee_i = squeeze(sumTerm); % msee of i-th stream
else % If H is n-dim, compute SINR values using H as page
    SSqrPageT = pagetranspose(SSqr); % page transpose this
    nVarVec = nVar*ones(1,size(W,2)); % create a vector to match SSqrPageT
    diagTerm = nVar./(SSqrPageT + nVarVec + eps); % eps to prevent divide by 0
    sumTerm = sum(absVSqr.*diagTerm,2);
    msee_i = permute(sumTerm,[3, 1, 2]); % msee of i-th stream
end
sinr = 1./msee_i - 1;
end