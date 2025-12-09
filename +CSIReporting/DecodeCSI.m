function [modulation,targetCodeRate,precoder] = DecodeCSI(carrier,pdsch,...
    pdschX,csiReport,csiFeedbackOpts)
%DecodeCSI decodes CSI to obtain new modulation, target code rate and
%   precoding matrix

csiReportConfig = csiFeedbackOpts.CSIReportConfig;
switch csiReportConfig.Mode
    case 'AUTOENCODER'
        % Using the decoder portion of the autoencoder to decode CSI
        [modulation,targetCodeRate,precoder] = decodeViaAutoencoder(...
            csiFeedbackOpts,csiReport,carrier,pdsch,pdschX);

    case 'RI-PMI-CQI'
        % Using 3GPP TS 38.211/214 RI-PMI-CQI decoding routine
        [modulation,targetCodeRate,precoder] = decodeViaRIPMICQI(...
            csiReportConfig,csiReport,carrier,pdschX);
end
end

%% Local Helper Fcns
function [modulation,targetCodeRate,precoder] = decodeViaRIPMICQI(...
    csiReportConfig,csiReport,carrier,pdschX)
%decodeRIPMICQI Subroutine to decode CSI using RI-PMI-CQI

% Map CSI to MCS
[modulation,targetCodeRate] = mapCSI2MCS(csiReportConfig.CQITable,csiReport);

% Map codebook-based precoding matrices from subbands to PRGs
precoder = mapPMISubband2PRGPrecodingMatrix(carrier,...
    pdschX.PRGBundleSize,csiReportConfig,csiReport.Precoder);
end

function [modulation,targetCodeRate,precoder] = decodeViaAutoencoder(...
    csiFeedbackOpts,csiReport,carrier,pdsch,pdschX)
%decodeAutoencoder Subroutine to decode CSI using the decoder portion of
%   the autoencoder

% Load the decoder body
decNet  = csiFeedbackOpts.CSIReportConfig.Autoencoder.Decoder;
aenOpts = csiFeedbackOpts.CSIReportConfig.Autoencoder.Options;

% Using the decoder body, do inference on the channel matrix as the
% codeword
H = Autoencoder.Decode(decNet,csiReport.H,aenOpts);

% Replicate channel matrix in the time domain to span 1
% slot worth of symbols
% Start by expanding a singleton dim in 2nd dim as placeholder and move 
% nRx,nTx to 3rd & 4th dim respectively
HTmp = permute(H,[1 4 2 3]); 
H = repelem(HTmp,1,carrier.SymbolsPerSlot,1,1,1); % repeat numSymbol along 2nd dim

% Compute PDSCH Tx parameters: modulation, target coderate and tx precoder
[modulation,targetCodeRate,precoder] = computePDSCHTxParameters(carrier,...
    pdsch,pdschX,H,csiReport.nVar);
end

function [modulation,tcr,wtx] = computePDSCHTxParameters(carrier,pdsch,pdschX,Hest,nVar)
%computePDSCHTxParameters compute MCS and MIMO precoder

% Get max rank based on Hest
maxRank = min(size(Hest,[3,4]));

% Allocate precoding matrix
WTmp = cell(maxRank,1);

% Allocate efficiency and mcsRank to NaN
eff = NaN(maxRank,1);
mcsRank = NaN(maxRank,2);

% Allocate MCS information structure
mcsInfo = struct();
mcsInfo.TableRow = [];
mcsInfo.TransportBLER = [];
mcsInfo = repmat(mcsInfo,1,maxRank);

% Calculate SINRs after precoding for each possible rank
for rnk = 1:maxRank
    % Perform SVD on the channel estimate
    w = computeSVDPrecodingMatrix(carrier,rnk,pdsch.PRBSet,Hest,nVar,pdschX.PRGBundleSize);

    % Store the MIMO precoding matrix w
    WTmp{rnk} = w;

    % Update pdsch
    pdsch.NumLayers = rnk;
    pdsch.ReservedRE = [];
    pdschX.XOverhead = 0;

    % Calculate SINR and select the corresponding best MCS
    [mcs,tmpInfo] = CSIReporting.SelectMCS(carrier,pdsch,pdschX,w,Hest,nVar);
    mcsInfo(rnk) = tmpInfo;

    % Get wideband MCS index and check if it's nonzero and nonnan
    mcsWb = mcs(1,:);
    nonzero = any(mcsWb);
    nonnan  = ~any(isnan(mcsWb));

    % Get number of codewords and fill in mcsRank array
    numCodewords = numel(mcsWb);
    mcsRank(rnk,:) = mcs(1,:);

    % Check if wideband mcs is nonzero and nonnan, then calculate
    % efficiency
    if nonzero && nonnan
        % Calculate throughput-related metrics using number of layers, code
        % rate, modulation and system-level BLER
        blerWb = mcsInfo(rnk).TransportBLER(1,:);
        numCodewordLayers = floor((rnk + (0:numCodewords-1))/numCodewords);
        
        % Update MCS
        mcs = mcsInfo(rnk).TableRow;

        % Calculate efficiency
        eff(rnk) = numCodewordLayers.*(1 - blerWb)*mcs(:,4);
    else
        eff(rnk) = 0;
    end
end

% Select the rank that maximizes throughput
[~,selectedRank] = max(eff);

% Update the PDSCH number of layers of the new rank
numLayers = selectedRank;

% Update MCS based on MCS index
numCodewords = ceil(numLayers/4); % 1 codeword can support up to 4 layers
mcs = mcsInfo(selectedRank).TableRow;

% Update modulation
Qm  = mcs(:,2).';
modLists = repmat({'QPSK','16QAM','64QAM','256QAM'}.',1,numCodewords);
modulation = modLists(Qm == repmat([2, 4, 6, 8].',1,numCodewords));

% Update target code rate
tcr = mcs(:,3).'/1024;

% Update Tx precoding matrix
wtx = WTmp{selectedRank};
end

function [modulation,tcr] = mapCSI2MCS(cqiTableNameInput,csiReport)

% Configure MCS based on CQI
cqi = csiReport.CQI(1,:); % Wideband
% Map RI to numCodewords according to TS 38.214 Table 5.1.3.1-1
% RI 1-4: 1 codeword
% RI 5-8: 2 codewords
numCodewords = ceil(csiReport.RI/4); 
% map CQI 0 to CQI 1 to prevent CQI from being 0. If CQI < 1, replace with
% 1, if CQI >= 1, keep as is
cqi = max([ones(1,numCodewords); cqi],[],1);

persistent thisTable cqiTableName
if isempty((thisTable)) || ~strcmpi(cqiTableName,cqiTableNameInput)
    cqiTableObj = nrCQITables;
    cqiTableName = cqiTableNameInput;
    thisTable = cqiTableObj.(['CQI',cqiTableName]);
end

modulation = thisTable.Modulation(cqi+1); % add 1 due to 0-based index
% Remove 'Out of Range' row
modulation = modulation(~strcmpi(modulation,'Out of Range'));
tcr  = thisTable.TargetCodeRate(cqi+1);
end

function wtx = mapPMISubband2PRGPrecodingMatrix(carrier,prgBundleSize,reportConfig,W)
%mapPMISubband2PRGPrecodingMatrix subroutine to map codebook-based 
%   precoding matrices from subbands to corresponding PRGs

if isempty(reportConfig.NSizeBWP)
    reportConfig.NSizeBWP = carrier.NSizeGrid;
end
if isempty(reportConfig.NStartBWP)
    reportConfig.NStartBWP = carrier.NStartGrid;
end

subbandInfo = CSIReporting.GetDLPMISubbandInfo(reportConfig);
wtx = CSIReporting.ComputePRGPrecoders(carrier,prgBundleSize,W,subbandInfo.SubbandSet);
wtx = permute(wtx,[2,1,3]);
end

function [wtx,sinr] = computeSVDPrecodingMatrix(carrier,numLayers,prbSet,HestGrid,nVar,prgBundleSize)
%computeSVDPrecodingMatrix calculates precoding matrices for all PRGs in 
% the carrier that overlap with the PDSCH allocation

persistent pdsch
if isempty(pdsch)
    pdsch = nrPDSCHConfig;
end

pdsch.NumLayers = numLayers;
pdsch.PRBSet = prbSet;
[wtx,hest] = CSIReporting.ComputeSVDPrecoders(carrier,pdsch,HestGrid,prgBundleSize);

nPRG = size(wtx,3); 
sinr = zeros(numLayers,nPRG);
for idx = 1:nPRG
    hestTmp = hest(:,:,idx);
    wtxTmp  = wtx(:,:,idx).';
    sinr(:,idx) = computePrecodedSINR2DHSubroutine(hestTmp,wtxTmp,nVar);
end
end

function sinr = computePrecodedSINR2DHSubroutine(H,W,nVar)
%computePrecodedSINR2DHSubroutine subroutine to compute precoded SINR with
%   channel matrix H being strictly 2D

% Reference
% nrPrecodedSINR(H,nVar,W) from Mathworks
% Li, Ping, et al. "On the Distribution of SINR for the MMSE MIMO Receiver 
% and Performance Analysis." IEEE Transactions on Information Theory, 
% vol. 52, no. 1, 1 Jan. 2006, pp. 271–286, 
% https://doi.org/10.1109/tit.2005.860466. Accessed 23 Sept. 2023.

assert(length(size(H)) == 2, 'computePrecodedSINR2DHSubroutine::H is not 2D matrix');

R = H*W;
[~,S,V] = pagesvd(R,'econ','vector');
SSqr = S.*S;
absVSqr = abs(V).^2;
diagTerm = nVar./(nVar + SSqr + eps); % eps to prevent divide by 0
diagTermPageT = pagetranspose(diagTerm);
sumTerm = sum(absVSqr.*diagTermPageT,2);
msee_i = squeeze(sumTerm); % msee of i-th stream
sinr = 1./msee_i - 1;
end