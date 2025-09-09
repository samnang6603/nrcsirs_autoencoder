function [CQI,PMISet,CQIInfo,PMIInfo] = CQISelect(carrier,CQIPMICalcParams,varargin)
%CQISelect Calculates CQI

[CQISubbandInfo,PMISubbandInfo] = getDLCSISubbandInfo(reportConfig);

% Get various parameters from CQIPMICalcParams
dmrsConfig = CQIPMICalcParams.DMRSConfig;
csirs = CQIPMICalcParams.CSIRS;
nTxAnts = unique(csirs.NumCSIRSPorts); % unique due to possibly more than one type of CSI-RS

% Put CDM type in a cell (possibly more than 1)
if iscell(csirs.CDMType)
    cdmType = csirs.CDMType;
else
    cdmType = {csirs.CDMType}; % if only 1
end

end

function [cqiSubbandInfo,pmiSubbandInfo] = getDLCSISubbandInfo(reportConfig)
% Returns downlink subband info on CQI and PMI
NSBPRB = reportConfig.SubbandSize;
reportConfig.CQISubbandSize = NSBPRB;
reportConfig.PMISubbandSize = NSBPRB;

% If have PRG size, then subband size is PRG size
% PRG: Precoding Resource block Group
if isempty(reportConfig.PRGSize)
    % If PRG isn't specified then do not ignore BWP
    reportConfig.IgnoreBWPSize = false;
else
    reportConfig.PMIMode = 'Subband';
    reportConfig.PMISubbandSize = reportConfig.PRGSize;
    % If PRG is specified, safely ignore BWP
    reportConfig.IgnoreBWPSize = true; 
end
cqiSubbandInfo = getSubbandInfo(reportConfig,'CQI'); 
pmiSubbandInfo = getSubbandInfo(reportConfig,'PMI');
end

function sbInfo = getSubbandInfo(reportConfig,Mode)
switch Mode
    case 'CQI'
        ignoreBWPSize = false; % Do not ignore BWP size for CQI
        NSBPRB = reportConfig.CQISubbandSize;
    case 'PMI'
        ignoreBWPSize = reportConfig.IgnoreBWPSize; % Dependent of PRG Size
        NSBPRB = reportConfig.PMISubbandSize;
    otherwise
        error('Unknown Reporting Mode')
end

nSizeBWP  = reportConfig.NSizeBWP;
nStartBWP = reportConfig.NStartBWP;

if strcmpi(reportingMode,'Wideband') || (~ignoreBWPSize && nSizeBWP < 24)
    % TS 38.214 Table 5.2.1.4-2: If nSizeBWP < 24 PRBs, then the number of
    % subbands is unity and the size is equal to nSizeBWP
    nSB = 1;
    sbSizes = nSizeBWP;
else
    % Compute first subband's size
    firstSBSize = NSBPRB - mod(nStartBWP,NSBPRB);

    % Compute final subband's size
    remdivNSBPRB = mod(nStartBWP + nSizeBWP,NSBPRB);
    if remdivNSBPRB ~= 0 % not divisible by NSBPRB
        finalSBSize = remdivNSBPRB; % then size is that remainder
    else
        finalSBSize = NSBPRB; % if divisible by the final size is NSBPRB
    end

    % Get total number of subbands, whcih is the total size of BWP,
    % excluding the first and final SB, divide by NSBPRB and then including
    % the first and final SB
    nSB = (nSizeBWP - firstSBSize - finalSBSize)/NSBPRB + 2;

    % Make a vector with each element representing the size of an SB
    sbSizes = NSBPRB*ones(1,nSB);
    sbSizes(1) = firstSBSize;
    sbSizes(end) = finalSBSize;
end
sbInfo.NumSubbands = nSB;
sbInfo.SubbandSizes = sbSizes;
end