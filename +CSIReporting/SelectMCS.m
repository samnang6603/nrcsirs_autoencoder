function [MCSIdx,MCSInfo] = SelectMCS(carrier,pdsch,pdschX,W,H,varargin)
%SelectMCS selects PDSCH modulation and coding scheme (MCS) based on
%   carrier, pdsch configuration, MIMO precoding matrx and estimated
%   channel response.

% Update PDSCH with number of layers from selected precoding matrix
pdsch.NumLayers = size(W,1);

% Get the start and size of the BWP (the grid)
NStartBWP = carrier.NStartGrid;
NSizeBWP  = carrier.NSizeGrid;
if ~isempty(pdsch.NStartBWP) || ~isempty(pdsch.NSizeBWP) % if BWP specified in PDSCH config
    NStartBWP = pdsch.NStartBWP;
    NSizeBWP  = pdsch.NSizeBWP;
end

% Noise Variance and BLER Target
if isempty(varargin)
    nVar = 1e-10; % default noise power
    targetBLER = 0.1;
elseif isscalar(varargin)
    nVar = varargin{1};
    targetBLER = 0.1;
else
    nVar = varargin{1};
    targetBLER = varargin{2};
end

% Get the start and end indices of the BWP and extract the H from that
% range
startIdx = 12*(NStartBWP - carrier.NStartGrid);
endIdx = startIdx + 12*NSizeBWP;
H = H(startIdx+1:endIdx,:,:,:);
[~,~,numRx,numTx] = size(H);

% Rearrange dim of precoding matrix, which is [NumLayers x Pcsirs x NPRG] 
% to [Pcsirs x NumLayers x NPRG] to help with SINR calculations
W = permute(W,[2, 1, 3]);

% If W has multiple precoding group, map each precoding matrix to its
% respective REs to calculate the SINR per RE
[numSubcarriers,numSyms,numPRGs] = size(H);
if numPRGs > 1
    W = mapMIMOPrecoder2RE(W,numSubcarriers,numSyms,...
        pdschX.PRGBundleSize,NStartBWP,NSizeBWP);
end

% Calculate SINR per RE
HRE = permute(H,[3,4,1,2]); % Rearrange dim to [nRx x nTx x nREs x nSyms]
HRE = reshape(HRE,numRx,numTx,[]); % Reshape to collapse numREs with numSyms
WRE = W(:,:,:);
SINRPerRE = computePrecodedSINRSubroutine(HRE,WRE,nVar);

% Lookup wideband MCS, effective SINR and BLER per subband
[mcsIdx,mcsTableRow,BLER] = L2SMapping.SelectMCSFromSystemLevel(carrier,...
    pdsch,pdschX.XOverhead,SINRPerRE,pdschX.MCSTable,targetBLER);

% Get outputs
MCSIdx = mcsIdx;
MCSInfo.TableRow = mcsTableRow;
MCSInfo.TransportBLER = BLER;

end

function WRE = mapMIMOPrecoder2RE(W,nSCs,nSyms,prgBundleSize,nStartBWP,nSizeBWP)
%mapMIMOPrecoder2RE Maps each precoding matrix page to its respective REs
% W is [Pcsirs x NumLayers x NPRG]

% Find the number of subbands and size of each subband
sbInfo = getSubbandInfo(nStartBWP,nSizeBWP,prgBundleSize);

% Allocate W per RE
WRE = zeros([size(W,1,2),nSCs,nSyms]);
sbStartIdx = 0;
for sbIdx = 1:sbInfo.NumSubbands
    thisSBSize = sbInfo.SubbandSizes(sbIdx);
    nRE = thisSBSize*12; % 12 subcarriers
    
    % RE indices
    reIdx = (1:nRE) + sbStartIdx;

    % Replicate the precoder matrix for all REs in the same subband
    thisSBPrecoder = W(:,:,sbIdx);
    repDims = [1,1,nRE,nSyms]; % repeat in this dim
    repArray = repmat(thisSBPrecoder,repDims);
    WRE(:,:,reIdx,:) = repArray;

    % Update the subband start index for next subband
    sbStartIdx = sbStartIdx + nRE;
end
end

function sbInfo = getSubbandInfo(nStartBWP,nSizeBWP,nSBPRB)
%getSubbandInfo Gets subband information of the bandwidth part

% Compute first subband's size
firstSBSize = nSBPRB - mod(nStartBWP,nSBPRB);

% Compute final subband's size
remdivNSBPRB = mod(nStartBWP + nSizeBWP,nSBPRB);
if remdivNSBPRB ~= 0 % not divisible by NSBPRB
    finalSBSize = remdivNSBPRB; % then size is that remainder
else
    finalSBSize = nSBPRB; % if divisible by the final size is NSBPRB
end

% Get total number of subbands, whcih is the total size of BWP,
% excluding the first and final SB, divide by NSBPRB and then including
% the first and final SB
numSB = (nSizeBWP - firstSBSize - finalSBSize)/nSBPRB + 2;

% Make a vector with each element representing the size of an SB
sbSizes = nSBPRB*ones(1,numSB);
sbSizes(1) = firstSBSize;
sbSizes(end) = finalSBSize;

% Output structure
sbInfo.NumSubbands = numSB;
sbInfo.SubbandSizes = sbSizes;
end

function sinr = computePrecodedSINRSubroutine(H,W,nVar)
%computePrecodedSINRSubroutine subroutine to compute precoded SINR

% Reference
% nrPrecodedSINR(H,nVar,W) from Mathworks
% Li, Ping, et al. "On the Distribution of SINR for the MMSE MIMO Receiver 
% and Performance Analysis." IEEE Transactions on Information Theory, 
% vol. 52, no. 1, 1 Jan. 2006, pp. 271–286, 
% https://doi.org/10.1109/tit.2005.860466. Accessed 23 Sept. 2023.

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
