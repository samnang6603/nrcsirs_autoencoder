function [PMISet,PMIInfo] = SelectDLPMI(carrier,csirs,csirsIndSubs,reportConfig,numLayers,H,nVar)
%DLPMISelect PMI Calculation for Downlink (PDSCH)

% Source of values:
%   - 3GPP TS 38.214 Table 5.2.2.2.1-2 (Type I single-panel codebook)
% Inspired by:
%   - MathWorks 5G Toolbox (hDLPMISelect.m), refactored into table form
%     for clarity, maintainability, and bug prevention.
%

% Configure Codebook
reportConfig = configureCodebook(reportConfig,csirs);

% Lookup codebook
[codebook,codebookIdxSetSizes] = lookupCodebook(reportConfig,csirs.NumCSIRSPorts,numLayers);

% Get PMI subband related information
subbandInfo = CSIReporting.GetDLPMISubbandInfo(reportConfig);

% Get CSI-RS indices inside BWP
[Hbwp,csirsIndBWP_k,csirsIndBWP_l,csirsIndBWP_p] = getInsideBWPComponents(carrier,...
    reportConfig,H,csirsIndSubs);

% Allocate NaN for PMI and get the information
[PMINaNSet,nanInfo] = allocatePMINaN(reportConfig,subbandInfo,codebook,...
    codebookIdxSetSizes,csirs.NumCSIRSPorts,numLayers,csirsIndBWP_k);

% If no CSI-RS Report all as NaN
% to be implemented later

% Compute SINRPerRE
Htmp = reshape(Hbwp,[],size(Hbwp,3),size(Hbwp,4));
Hcsirs = Htmp(csirsIndBWP_k+(csirsIndBWP_l-1)*size(Hbwp,1),:,:);
SINRPerRE = CSIReporting.ComputePrecodedSINRPerRE(Hcsirs,codebook,...
    codebookIdxSetSizes,nVar,csirsIndBWP_k,numLayers);

% Compute PMI Wideband
[PMISet,W,SINRPerREPMI,subbandSINRs] = Codebook.ComputeTypeIPMIWideband(SINRPerRE,...
    codebook,codebookIdxSetSizes);

% Compute PMI Subband
numSubbands = subbandInfo.NumSubbands;
if numSubbands > 1
    switch reportConfig.CodebookType
        case {'TypeISinglePanel','TypeIMultiPanel'}
            [PMISet,W,SINRPerREPMI,subbandSINRs] = Codebook.ComputeTypeIPMISubband(...
                reportConfig,codebook,PMISet,SINRPerRE,subbandInfo,...
                codebookIdxSetSizes,numLayers,csirsIndBWP_k);
        case {'TypeII','eTypeII'}
    end
end

% Output result
switch reportConfig.CodebookType
    case {'TypeISinglePanel','TypeIMultiPanel'}
        SINRPerRE_out = SINRPerRE; % for all layers for all PMI indices
        codebook_out  = codebook; % codebook for all PMI indices
        % SINR per subband for all layers and all PMI indices
        subbandSINRs = reshape(subbandSINRs,[subbandInfo.NumSubbands,...
            numLayers,codebookIdxSetSizes]);
    case {'TypeII','eTypeII'}
        SINRPerRE_out = [];
        codebook_out  = [];
end
PMIInfo.SINRPerRE = SINRPerRE_out;
PMIInfo.SINRPerREPMI = SINRPerREPMI;   % SINR value per RE for all the layers for the reported PMI
PMIInfo.Codebook = codebook_out;
PMIInfo.SINRPerSubband = subbandSINRs; % SINR value per subband for all the layers for the reported PMI
PMIInfo.W = W;                         % Precoding matrix corresponding to the reported PMI
PMIInfo.CSIRSIndices = [csirsIndBWP_k, csirsIndBWP_l, csirsIndBWP_p]; % CSI-RS RE indices where the SINR value per RE are calculated

end

function reportConfig = configureCodebook(reportConfig,csirs)
%configureCodebook Configure codebook and validate panel dimensions

% Start with N1, N2, O1, O2 = 1
numCSIRSPorts = csirs.NumCSIRSPorts(1);

% TypeISinglePanelFlag = strcmpi(reportConfig.CodebookType,'TypeISinglePanel');
% TypeIMultiPanelFlag  = strcmpi(reportConfig.CodebookType,'TypeIMultiPanel');
% TypeIIFlag  = strcmpi(reportConfig.CodebookType,'TypeII');
% ETypeIIFlag = strcmpi(reportConfig.CodebookType,'ETypeII');

switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        reportConfig = Codebook.TypeISinglePanel.ConfigureCodebookParameters(reportConfig,numCSIRSPorts);
    case 'TypeIMultiPanel'
        % To be implemented later when moving to multiPanel
end





end


function [codebook,codebookIdxSetSizes] = lookupCodebook(reportConfig,numCSIRSPorts,numLayers)
%lookupCodebook Extract codebook. Calls codebook type-specific subroutine
%to configure, calculate and extract codebook.
%

switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        if numCSIRSPorts == 1
            % If the number of port is 1, then codebook is a scalar value
            % of unity
            codebook = 1;
        else
            % Codebook is an N-dim of size
            % Pcsirs x numLayers x i2Len x i11Len x i12Len x i13Len
            % - OR -
            % Pcsirs x numLayers x i2Len
            [codebook,codebookIdxSetSizes] = Codebook.TypeISinglePanel.ExtractCodebook(reportConfig,numCSIRSPorts,numLayers);
        end

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
bwpStartPos = reportConfig.NStartBWP - carrier.NStartGrid;

% Calculate the SINR and CQI values
% csirsInd = [subcarrier, symbol, port]
kTmp = csirsIndSubs(:,1); % temp subcarrier indices
lTmp = csirsIndSubs(:,2); % temp symbol indices

% Only consider the CSI-RS indices inside the BWP
% Since BWP is given in RBs, it has to be multiplied by NSC = 12 to
% convert the BWP to effective sc units.
% The CSI-RS to the right of the BWP is simply the region starting from the
% BWP SC to the right
csirsRightOfBWPStart = kTmp >= bwpStartPos*12 + 1;

% While the CSI-RS to the left of the BWP is simply the region starting
% from the end of the BWPStart + BWPSize to the left
csirsLeftOfBWPEnd  = kTmp <= (bwpStartPos + reportConfig.NSizeBWP)*12;

% Therefore, the intersection is the area where valid SCs containing CSI-RS
% exist inside BWP
indInsideBWP = csirsRightOfBWPStart & csirsLeftOfBWPEnd;

% Update the subcarrier and symbol indices
% Also convert the CSI-RS SCs subscripts to BWP scale
k = kTmp(indInsideBWP) - bwpStartPos*12;
l = lTmp(indInsideBWP);
p = ones(size(k)); % only the first CSI-RS port matters

% Channel estimate inside the BWP
% Only change the k index range w.r.t carrier
Hbwp = H(bwpStartPos*12 + 1:(bwpStartPos + reportConfig.NSizeBWP)*12,:,:,:);

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
    nanInfo.W = NaN(numCSIRSPorts,numLayers);
end

end
