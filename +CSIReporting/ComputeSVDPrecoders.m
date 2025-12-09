function [Wprg,Hprg] = ComputeSVDPrecoders(carrier,pdsch,HestGrid,prgSize)
%ComputeSVDPrecoder uses SVD decomposition to compute precoding matrix from
%   channel estimation Hest
%   Hest is shaped (12*NSizeGrid x nSym x nRx x nTx)

% Get Precording Resource Group info
prgInfo = nrPRGInfo(carrier,prgSize);

% Set the grid for RB relative to NStartGrid
gridRBSet = setGridRBSet(carrier,pdsch);

% Get the dimensions of Hest
[K,nSym,nRx,nTx] = size(HestGrid);

% Reshape the the first dim of H (12sc*NSizeGrid) to (12sc x NsizeGrid)
HestGrid = reshape(HestGrid,[12, K/12, nSym, nRx, nTx]);

% Take the average of the channel estimate across all subcarriers (1st dim)
% and then across all OFDM symbols (3rd dim)
HestAvgSub = mean(HestGrid,1);
HestAvgSym = mean(HestAvgSub,3);

% Extract all allocated RBs
% HestRB is 1(AvgSub) x gridRBSet x 1(AvgSym) x nRx x nTx
HestRB = HestAvgSym(:,gridRBSet + 1,:,:,:);

% Permute to shape nRx x nTx x Allocated_PRBs
HestRB = squeeze(HestRB); % Squeeze singleton dims (1st and 3rd)
HestRB = permute(HestRB,[2, 3, 1]);

% Per PRG
numLayers = pdsch.NumLayers;
Wprg = zeros(numLayers,nTx,prgInfo.NPRG);
Hprg = zeros(nRx,nTx,prgInfo.NPRG);
pdschPRGs = prgInfo.PRGSet(gridRBSet + 1);

for idx = unique(pdschPRGs).'
    % Take the average of the channel estimate across all allocated RBs
    % except zero estimates for incomplete PRGs or missing information
    Hprb = HestRB(:,:,pdschPRGs == idx);
    nonzeroIdx = ~all(Hprb == 0,[1, 2]);
    HprgTmp = mean(Hprb(:,:,nonzeroIdx),3); % 3rd dim is allocated RBs for PRGs
    Hprg(:,:,idx) = HprgTmp;

    % Do SVD and extract only the right singular vectors matrix V due to
    % transmission side
    [~,~,V] = svd(HprgTmp); % right singular vectors matrix
    WTmp = permute(V(:,1:numLayers),[2, 1, 3]);
    Wprg(:,:,idx) = WTmp/sqrt(numLayers); % normalize
end
end


function gridRBSet = setGridRBSet(carrier,pdsch)
%arrangeGridRBSet get and arrange allocated RBs in the carrier resource
%   grid relative to the start of the grid, NStartGrid

if (pdsch.VRBToPRBInterleaving)
    [~,indinfo] = nrPDSCHIndices(carrier,pdsch);
    gridRBSet = indinfo.PRBSet;
else
    gridRBSet = pdsch.PRBSet;
end

if isempty(pdsch.NStartBWP)
    % If empty, start fresh from 0
    bwpOffset = 0;
else
    % If NStartBWP is specified then find offset relative to it
    bwpOffset = pdsch.NStartBWP - carrier.NStartGrid;
end
gridRBSet = gridRBSet + bwpOffset;

end