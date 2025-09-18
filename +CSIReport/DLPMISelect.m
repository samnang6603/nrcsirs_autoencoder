function [PMISet,PMIInfo] = DLPMISelect(carrier,csirs,csirsIndSub,reportConfig,numLayers,H,nVar)
%DLPMISelect PMI Calculation for Downlink (PDSCH)

% Source of values:
%   - 3GPP TS 38.214 Table 5.2.2.2.1-2 (Type I single-panel codebook)
% Inspired by:
%   - MathWorks 5G Toolbox (hDLPMISelect.m), refactored into table form
%     for clarity, maintainability, and bug prevention.
%

reportConfig = Codebook.TypeISinglePanel.ConfigureCodebookParameters(reportConfig);

% Lookup codebook
[codebook,codebookIdxSetSizes] = getCodebook(reportConfig,csirs.NumCSIRSPorts,numLayers);

% Get PMI subband related information
subbandInfo = CSIReport.DLPMISubbandInfo(carrier,reportConfig);

% Get CSI-RS indices inside BWP
[H_bwp,csirsIndBWP_k,csirsIndBWP_l,csirsIndBWP_p] = getInsideBWPComponents(carrier,reportConfig,H,csirsIndSubs);

% Allocate NaN for PMI and get the information
[PMINaNSet,nanInfo] = allocatePMINaN(reportConfig,subbandInfo,codebook,codebookIdxSetSizes,csirs.numCSIRSPorts,numLayers,csirsIndBWP_k);

% Compute PMI
Htmp = reshape(Hbwp,[],size(Hbwp,3),size(Hbwp,4));
Hcsirs = Htmp(csirsIndBWP_k+(csirsIndBWP_l-1)*size(H_bwp,3),:,:);
SINRPerRE = CSIReport.ComputePrecodedSINRPerRE(Hcsirs,codebook,...
    codebookIdxSetSizes,nVar,csirsIndBWP_k,numLayers);




end



function [codebook,codebookIdxSetSizes] = getCodebook(reportConfig,numCSIRSPorts,numLayers)
%
%

switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        if numCSIRSPorts == 1
            % If the number of port is 1, then codebook is a scalar value
            % of unity
            codebook = 1;
        else
            % Codebook is an N-dim of size
            % nHEstCSIRSPorts x numLayers x i2Len x i11Len x i12Len x i13Len
            % - OR -
            % nHEstCSIRSPorts x numLayers x i2Len
            reportConfig = Codebook.TypeISinglePanel.ConfigureCodebookParameters();
            [codebook,codebookSize] = Codebook.TypeISinglePanel.ExtractCodebook(reportConfig,numCSIRSPorts,numLayers);
        end

        % Size of the codebook
        i2Len  = codebookSize(3);
        i11Len = codebookSize(4);
        i12Len = codebookSize(5);
        i13Len = codebookSize(6);

        codebookIdxSetSizes = [i2Len; i11Len; i12Len; i13Len];

    case 'TypeIMultiPanel'

    case 'TypeII'

end
end

function [Hbwp,k,l,p] = getInsideBWPComponents(carrier,reportConfig,H,csirsIndSubs)
%getInBWPComponents
%
% Within a PRB, find:
% k: subcarrier index (frequency axis)
% l: OFDM symbol index (time axis)
% p: CSI-RS port index (antenna port)
% H: Channel estimated inside the BWP (K-L-nRxAnts-Pcsirs)
%
% Converts CSI-RS indices and channel estimates from
% carrier-level coordinates to BWP-relative coordinates. It extracts
% only the REs and corresponding channel estimates that fall within
% the UE's active BWP, allowing downstream CQI/PMI/SINR calculations
% to be computed accurately over the allocated bandwidth.
%
% The output indices (k,l,p) can be used to directly index into Hbwp
% to access the per-port channel estimate for each CSI-RS RE.
%
% Why H (channel estimate) needs to be restricted to the BWP:
%   Hbwp contains only the subcarriers of interest for the UE:
%     - Ensures CQI/PMI are computed only on allocated subcarriers
%     - Prevents inclusion of irrelevant carrier subcarriers which would
%       skew SINR calculations
%     - Reduces computational load by slicing out only necessary
%       subcarriers

% Find the start of BWP relative to the carrrier
bwpStart = reportConfig.NStartBWP - carrier.NStartGrid;

% Calculate the SINR and CQI values
% csirsInd = [subcarrier, symbol, port]
kTmp = csirsIndSubs(:,1); % temp subcarrier indices
lTmp = csirsIndSubs(:,2); % temp symbol indices

% Only consider the CSI-RS indices inside the BWP
% Since BWP is given in RBs, it has to be multiplied by NSC = 12 to
% convert the BWP to effective sc units.
% The CSI-RS to the right of the BWP is simply the region starting from the
% BWP SC to the right
csirsRightOfBWPStart = kTmp >= bwpStart*12 + 1;

% While the CSI-RS to the left of the BWP is simply the region starting
% from the end of the BWPStart + BWPSize to the left
csirsLeftOfBWPEnd  = kTmp <= (bwpStart + reportConfig.NSizeBWP)*12;

% Therefore, the intersection is the area where valid SCs containing CSI-RS
% exist inside BWP
indInsideBWP = csirsRightOfBWPStart & csirsLeftOfBWPEnd;

% Update the subcarrier and symbol indices
% Also convert the CSI-RS SCs subscripts to BWP scale
k = kTmp(indInsideBWP) - bwpStart*12;
l = lTmp(indInsideBWP);
p = ones(size(k)); % only the first CSI-RS port matters

% Channel estimate inside the BWP
% Only change the k index range w.r.t carrier
Hbwp = H(bwpStart*12 + 1:(bwpStart + reportConfig.NSizeBWP)*12,:,:,:);

end

function [PMINaNSet,nanInfo] = allocatePMINaN(reportConfig,subbandInfo,...
    codebook,codebookIdxSetSizes,numCSIRSPorts,numLayers,csirsIndBWP_k)

csirsIndLen = length(csirsIndBWP_k);
nanInfo.CSSIRSIndices = NaN(csirsIndLen,3);

switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        PMINaNSet.i1 = NaN(1,3);
        PMINaNSet.i2 = NaN(1,subbandInfo.NumSubbands);
    case 'TypeIMultiPanel'
        PMINaNSet.i1 = NaN(1,6);
        PMINaNSet.i2 = NaN(3,subbandInfo.NumSubbands);
end

if contains(reportConfig.CodebookType,'TypeI')
    nanInfo.SINRPerRE = NaN([csirsIndLen,numLayers,codebookIdxSetSizes]);
    nanInfo.SINRPerREPMI = NaN(csirsIndLen,numLayers);
    nanInfo.SINRPerSubband = NaN(subbandInfo.NumSubbands,numLayers);
    nanInfo.Codebook = codebook;
    nanInfo.W = NaN(numCSIRSPorts,nLayers);
end

end
