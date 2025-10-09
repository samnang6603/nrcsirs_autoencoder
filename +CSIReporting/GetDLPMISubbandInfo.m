function PMISubbandInfo = GetDLPMISubbandInfo(reportConfig)
%DLPMISubbandInfo Get information on downlink PMI subband
%
% This implementation is adapted from MathWorks example helper 
% hDLPMISubbandInfo.m. Original at MathWorks, re-written here 
% for research/educational use.

% Reference: 
% 1. hDLPMISubbandInfo() function as part of Mathwork's example "PDSCH
%    Throughput Using CSI Feedback".
%    URL: https://www.mathworks.com/help/5g/ug/nr-pdsch-throughput-using-csi-feedback.html

% Get the start of size of BandWidth Part
nStartBWP = reportConfig.NStartBWP; 
nSizeBWP  = reportConfig.NSizeBWP;

% Number of Subband per PRB
nSBPRB = reportConfig.SubbandSize;

if strcmpi(reportConfig.PMIMode,'Wideband') || (nSizeBWP < 24)
    % TS 38.214 Table 5.2.1.4-2 states that if BWP is less than 24 PRBs,
    % then we cannot divide BWP into subbands. The number of subbands is
    % unity in this case and the subband size = BWP size
    nSubbands = 1;
    nSBPRB = nSizeBWP;
    subbandSizes = nSBPRB;
    subbandSet = ones(1,nSizeBWP);
else
    % Divide BWP into subbands
    prb = nStartBWP + (0:nSizeBWP-1);
    prbInt = floor(prb/nSBPRB);
    subbandSet = 1 + prbInt;
    subbandSizes = histcounts(prbInt,BinMethod='integers');
    nSubbands = length(subbandSizes);
end

PMISubbandInfo.SubbandSizes = subbandSizes;
PMISubbandInfo.SubbandSet = subbandSet;
PMISubbandInfo.NumSubbands = nSubbands;

end