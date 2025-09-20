function [CQI,PMISet,CQIInfo,PMIInfo] = SelectCQI(CQIPMICalcParams)
%CQISelect Calculates CQI Based on the given CQI-PMI-RI configurations 

% Get various parameters from CQIPMICalcParams

% Carrier config
carrier = CQIPMICalcParams.Carrier;

% DMRS Config
dmrsConfig = CQIPMICalcParams.DMRSConfig;

% CSI-RS
csirs = CQIPMICalcParams.CSIRS;

% Get CSI-RS Indices and nTxAnts
csirsInd = getCSIRSIndices(carrier,csirs); % [sc, symbol, port]
numTxAnts = unique(csirs.NumCSIRSPorts); % unique due to possibly more than one type of CSI-RS

% Report Config
reportConfig = CQIPMICalcParams.ReportConfig;

% Allocate empty PRG Size if no PRG Size
reportConfig.PRGSize = [];

% Number of layer
% This is directly tied to the rank indicated by the UE (RI). Since we sift
% through each possible valid rank, the nLayer is defined by the current
% rank we're on.
numLayers = CQIPMICalcParams.ThisRank;

% Get the channel matrix
H = CQIPMICalcParams.Channel;

% Get the Noise variance
nVar = CQIPMICalcParams.NoiseVariance;

% Allocate empty SINRTable
SINRTable = [];
inputSINRTable = 0;

% Find the number of subbands and the size of each subband for the given
% CQI and PMI configuration.
[CQISubbandInfo,PMISubbandInfo] = getDLCSISubbandInfo(reportConfig);
CQINumSubbands = CQISubbandInfo.NumSubband;

% Number of codewords
% If nLayer > 4 -> 2 codewords
numCodewords = ceil(numLayers/4);

% Find the start of BWP relative to the carrrier
bwpStart = reportConfig.NStartBWP - carrier.NStartGrid;

% Calculate the SINR and CQI values
% csirsInd = [subcarrier, symbol, port]
csirsIndSubs_kTmp = csirsInd(:,1); % temp subcarrier indices
csirsIndSubs_lTmp = csirsInd(:,2); % temp symbol indices

% Only consider the CSI-RS indices inside the BWP
% Since BWP is given in RBs, it has to be multiplied by NSC = 12 to
% convert the BWP to effective sc units.
% The CSI-RS to the right of the BWP is simply the region starting from the
% BWP SC to the right
csirsRightOfBWPStart = csirsIndSubs_kTmp >= bwpStart*12 + 1;

% While the CSI-RS to the left of the BWP is simply the region starting
% from the end of the BWPStart + BWPSize to the left
csirsLeftOfBWPEnd  = csirsIndSubs_kTmp <= (bwpStart + reportConfig.NSizeBWP)*12;

% Therefore, the intersection is the area where valid SCs containing CSI-RS
% exist inside BWP
indInsideBWP = csirsRightOfBWPStart & csirsLeftOfBWPEnd;

% Update the subcarrier and symbol indices
% Also convert the CSI-RS SCs subscripts to BWP scale (-bwpStart*12)
csirsIndSubs_k = csirsIndSubs_kTmp(indInsideBWP) - bwpStart*12;
csirsIndSubs_l = csirsIndSubs_lTmp(indInsideBWP);

% PMI select using DLPMISelect()
[PMISet,PMIInfo] = CSIReporting.SelectDLPMI(carrier,csirs,csirsInd,reportConfig,numLayers,H,nVar);

% Handles case where there is no CSI-RS
% >>>>>>>>>>>>>>>>>> TBE

sinrPerREPMI = PMIInfo.SINRPerREPMI;
switch reportConfig.CodebookType
    case {'TypeISinglePanel','TypeIMultiPanel'}

        % Allocate SINRPerSubband with NaN
        SINRPerSubband = NaN(CQISubbandInfo.NumSubbands,numLayers);

        % PRGSize Handling
        if isfield(reportConfig,'PRGSize') && ~isempty(reportConfig.PRGSize)
            % PRG Bundle Size TBE
        else
            % IF PRGSize is not specified, the PMI selection is either in
            % wideband or subband level granularity

            % Compute SINR values of the PMISet in RB level granularity.
            % This is used for information purpose only.
            % >>>>>>>>>>>>>>>>>> TBE getSINRPerRB fcn
            % SINRPerRBPerCodeword
            
            switch reportConfig.PMIMode
                case 'Wideband'
                    % When PMI mode is set to Wideband granularity, only i2
                    % index is reported. The SINR values, SINRPerSubband,
                    % are calculated for the entire BWP span. 
                    SINRPerSubband = getSubbandSINR(sinrPerREPMI,CQINumSubbands,csirsIndSubs_k);
                case 'Subband'
                    % When PMI mode is set to Subband granularity, for
                    % codebook type:
                    % - TypeISinglePanel: one i2 value is reported per
                    %   subband.
                    % - TypeIMultiPanel: a set of [i20 i21 i22] are
                    %   reported per subband.
                    TypeISinglePanelFlag = strcmpi(reportConfig.CodebookType,'TypeISinglePanel');
                    for idx = 1:size(PMISet.i2,2)
                        if ~any(isnan(PMISet.i2(:,idx)))
                            if TypeISinglePanelFlag
                                tmp = PMIInfo.SINRPerSubband(idx,:,PMISet.i2(idx),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3));
                            else
                                tmp = PMIInfo.SINRPerSubband(idx,:,PMISet.i2(1,idx),PMISet.i2(2,idx),PMISet(3,idx),PMISet.i1(1),PMISet.i1(2),PMISet.i1(3),PMISet.i1(4),PMISet.i1(5),PMISet.i1(6));
                            end
                        end
                        SINRPerSubband(idx,:) = tmp;
                    end
            end
        end

    case {'TypeII','ETypeII'}
        % TBE
end

% Allocate SINR per subband per Codeword
SINRPerSubbandPerCodeword = zeros(CQINumSubbands,numCodewords);
for idx = 1:CQINumSubbands
    % Get the SINR values per layer
    SINRPerLayer = squeeze(SINRPerSubband(idx,:));
    
    % Calculate the SINR values per codeword
    if ~any(isnan(SINRPerLayer))
        layerDemapped = nrLayerDemap(SINRPerLayer); % demap the layer
        SINRPerCodeword = zeros(1,length(layerDemapped));
        for b = 1:length(layerDemapped)
            SINRPerCodeword(b) = sum(layerDemapped); % sum the SINR per layer
        end
    else
        % If there are NaN values, that means no CSI-RS resources present
        % in the current subband
        SINRPerCodeword = NaN(1,numCodewords);
    end
    SINRPerSubbandPerCodeword(idx,:) = SINRPerCodeword;
end




end


function csirsInd = getCSIRSIndices(carrier,csirsConfig)
% Get CSI-RS indices
% Put CDM type in a cell (possibly more than 1)
if iscell(csirsConfig.CDMType)
    cdmType = csirsConfig.CDMType;
else
    cdmType = {csirsConfig.CDMType}; % if only 1
end

if ~iscell(csirsConfig.CSIRSType)
    csirsConfig.CSIRSType = {csirsConfig.CSIRSType};
end

% Ignore ZP CSI-RS because it is not used for CSI estimation
% Find the number of resources occupied by ZP
numZP = sum(strcmpi(csirsConfig.CSIRSType,'zp'));

% Get CSI-RS indices
tmpInd = nrCSIRSIndices(carrier,csirsConfig,IndexStyle='subscript',...
    OutputResourceFormat='cell');
% Chop off the ZP indices becausse nrCSIRSIndices return the ZP indices 
% first, if any.
tmpInd = tmpInd(numZP+1:end)';

% Extract the first port's NZP CSI-RS indices
for nzpIdx = 1:length(tmpInd) % number of NZP type
    nzpInd = tmpInd{nzpIdx};
    port1Ind = nzpInd(:,3) == 1; % first port indices only
    tmpInd{nzpIdx} = nzpInd(port1Ind,:);
end

% Extract the indices with the lowest RE of each CDM group. This will limit
% the number of CSI-RS REs and speed up computation
if ~strcmpi(cdmType{1},'noCDM')
    for resourceIdx = 1:numel(tmpInd)
        totalIndices = size(tmpInd{resourceIdx},1);
        switch cdmType{1}
            case 'FD-CDM2'
                indPerSym = totalIndices;
            case 'CDM4'
                indPerSym = totalIndices/2;
            case 'CDM8'
                indPerSym = totalIndices/4;
            otherwise
                error('Invalid CDM type')
        end
        tmpIndPerSym = tmpInd{resourceIdx}(1:indPerSym,:);
        
        % Downsample by 2 regardless of CDM type because we only need one 
        % RE per CDM group to represent the group and to improve 
        % computational speed 
        tmpInd{resourceIdx} = tmpIndPerSym(1:2:end,:);
    end
end

if ~isempty(tmpInd)
    csirsInd = cell2mat(tmpInd);
else
    csirsInd = [];
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
        reportingMode = reportConfig.CQIMode;
    case 'PMI'
        ignoreBWPSize = reportConfig.IgnoreBWPSize; % Dependent of PRG Size
        NSBPRB = reportConfig.PMISubbandSize;
        reportingMode = reportConfig.PMIMode;
    otherwise
        error('Unknown Reporting Mode')
end

nSizeBWP  = reportConfig.NSizeBWP;
nStartBWP = reportConfig.NStartBWP;

% Valid BWP status if Size < 24 and is considered
validBWPStatus = ~ignoreBWPSize && nSizeBWP < 24;

if strcmpi(reportingMode,'Wideband') || validBWPStatus
    % TS 38.214 Table 5.2.1.4-2: If nSizeBWP < 24 PRBs, then the number of
    % subbands is unity and the size is equal to nSizeBWP
    numSB = 1;
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
    numSB = (nSizeBWP - firstSBSize - finalSBSize)/NSBPRB + 2;

    % Make a vector with each element representing the size of an SB
    sbSizes = NSBPRB*ones(1,numSB);
    sbSizes(1) = firstSBSize;
    sbSizes(end) = finalSBSize;
end
sbInfo.NumSubbands = numSB;
sbInfo.SubbandSizes = sbSizes;
end

function SINRPerSubband = getSubbandSINR(SINRPerREPMI,NumSubbands,k)
%getSubbandSINR Calculate the average SINR values throughout all REs within
%the subband span of one slot, with respect to the reported PMI indices.

% Allocate using NaN
subbandSINRs = NaN(NumSubbands,size(SINRPerREPMI,2));
subbandStartPos = 0;
for thisSB = 1:NumSubbands
    subbandSize = SubbandSizes(thisSB);
    % Lower bound of BWP
    lowerBound = k >= (subbandStartPos*12 + 1);

    % Upper bound requires sb idx
    upperBound = k <= ((subbandStartPos + subbandSize)*12);

    % Intersection between lower and upper bound
    subbandInd = lowerBound & upperBound;

    % SINR
    SINR = SINRPerREPMI(subbandInd,:,:);
    if ~all(isnan(SINR(:)))
        avgSINR = mean(SINR,1);
        subbandSINRs(thisSB,:) = avgSINR;
    end
    % Update subband start position
    subbandStartPos = subbandStartPos + subbandSize;
end
end

function SINRPerRBPerCodeword = getSINRPerRB()


end
 