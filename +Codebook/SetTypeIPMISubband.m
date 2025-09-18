function [PMISet,W,SINRPerREPMI,subbandSINRs] = SetTypeIPMISubband(carrier,...
    codebook,PMISet,SINRPerRE,subbandInfo,codebookIdxSetSizes,numLayers,k,l)
%SetTypeIPMISubband Computes PMI for Type I single or multipanel for
%subband. Takes results from wideband computations as well

numSubbands = subbandInfo.NumSubbands;
numCSIRSPorts = size(codebook,1);

% Allocations some outputs
W = zeros(numCSIRSPorts,numLayers,numSubbands);
SINRPerREPMI = zeros(length(k),numLayers);
subbandSINRs = NaN(numSubbands,numLayers,codebookIdxSetSizes);

% Start of subband is index 0
subbandInit = 0;

% Lower bound of BWP
lowerBound = k >= (subbandInit*12 + 1);


% Sift through the subbands
for idx = 1:numSubbands

    % Subband size w.r.t subband index
    subbandSize = subbandInfo.SubbandSizes(idx);
    
    % Upper bound requires sb idx
    upperBound = k <= ((subbandInit + subbandSize)*12 + 1);

    % Intersection between lower and upper bound
    subbandInd = lowerBound & upperBound;

    % SINR value per subband
    sinrValPerSubband = SINRPerRE(subbandInd,:,:,:,:,:,:,:,:,:,:);

end

end