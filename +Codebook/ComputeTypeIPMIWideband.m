function [PMISet,W,SINRPerREPMI,subbandSINRs] = ComputeTypeIPMIWideband(SINRPerRE,codebook,codebookIdxSetSizes)
%SetPMIWideband output PMI set for Type I single panel and multiplanel
%codebooks, precoding matrix and subband SINR info.

% Sum the SINR along BWP and then along all layer and collapse singleton
% dimensions
% SINRPerRE is size sc(BWP) x 1 x nLayer
SINRSum = sum(SINRPerRE,1); % sum along BWP
SINRSum = squeeze(sum(SINRSum,2)); % sum along layers and collapse singleton dim

% Reshape the SINR into codebookIdxSetSizes dim
SINRSum = reshape(SINRSum,codebookIdxSetSizes);

% Round SINR value to 4 decimal places to avoid fluctuations in the PMI
% output due to minute variations of the SINR values and their respective
% PMI indices
SINRSum = round(SINRSum,4,'decimals');

% Get the indices set that correspond to the precoding matrix with max SINR
maxSINR = max(SINRSum,[],'all');
maxSINRInd = find(totalSINR == maxSINR,1); % find maxSINR indices 1st dim
[i2, i11, i12, i13] = ind2sub(size(SINRSum),maxSINRInd);
PMISet.i1 = [i11, i12, i13];
PMISet.i2 = i2;
W = codebook(:,:,i2,i11,i12,i13); % Lookup and extract precoding matrix from PMISet
SINRPerREPMI = SINRPerRE(:,:,i2,i11,i12,i13); % Get the corresponding SINR from PMISet

% Get the average for Subband SINR
subbandSINRs = mean(SINRPerRE,1); 
end