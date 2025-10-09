function [Wprg,prgBands] = ComputePRGPrecoders(carrier,prgSize,WperBand,bandset)
%ComputePRGPrecoders Compute PRG precoders from precoders for arbitrary 
% carrier RBs

% Convert band-level precoders (WperBand) into PRG-level precoders (Wprg)
% Each PRG spans multiple RBs in frequency. For each PRG:
%   - Identify which frequency "band" it overlaps via bandset
%   - Select the lowest valid band index if multiple overlap
%   - Leave as NaN if the PRG has no defined band
%
% Once mapped, copy corresponding precoders from Wband into Wprg.
% Undefined PRGs (no valid band) are filled with zeros.
%
% Result:
%   Wprg(:,:,n) = precoder matrix for the nth PRG
%   Wprg aligns in frequency with the carrier's PRG structure.


% prgSize: number of PRG per RB
prgInfo = nrPRGInfo(carrier,prgSize);
prgSet  = prgInfo.PRGSet;
uPRG = unique(prgSet);
prgBands = zeros(length(uPRG),1);
for b = 1:length(uPRG)
    isThisPRGSet = prgSet == uPRG(b);
    prgBands(b) = min(bandset(isThisPRGSet));
end

% Concatenate the NaN array to the start of prgBands for any complete PRGs
% that don't start at the beginning of the carrier grid (maybe offset).
% This is necessary to line up 'Wprg' correctly when indexing later
nanArr = NaN(prgSet(1)-1,1);
prgBands = [nanArr;
            prgBands];

% Now create the final PRG-level precoder cube 'Wrpg'. Return zeros for
% undefined PRGs
Wprg = zeros([size(WperBand,1:2), prgInfo.NPRG]);
nonnanIdx = ~isnan(prgBands);
Wprg(:,:,nonnanIdx) = WperBand(:,:,prgBands(nonnanIdx));
end