function [Wprg,Hprg] = ComputeSVDPrecoder(carrier,pdsch,Hest,prgSize)
%ComputeSVDPrecoder uses SVD decomposition to compute precoding matrix from
% channel estimation Hest

% Get Precording Resource Group info
prgInfo = nrPRGInfo(carrier,prgSize);

% Set the grid for RB relative to NStartGrid
gridRBSet = setGridRBSet(carrier,pdsch);

% Get the dimensions of Hest
[K,L,R,P] = size(Hest);


end


function gridRBSet = setGridRBSet(carrier,pdsch)
%arrangeGridRBSet get and arrange allocated RBs in the carrier resource
% grid relative to the start of the grid, NStartGrid

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