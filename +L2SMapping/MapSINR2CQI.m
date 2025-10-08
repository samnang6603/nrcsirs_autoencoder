function [L2SMConfig,cqiIdx,cqiInfo] = MapSINR2CQI(L2SMConfig,carrier,phySigRes,...
    xOverhead,sinr,Qm,TCR,blerThresh)

% Check for new CQI Table, 
selectedCQIQmTCR = [Qm, TCR];
L2SMCQITable = L2SMConfig.CQI.Table;
isNewTable = isempty(L2SMCQITable) || ~isequaln(L2SMCQITable,selectedCQIQmTCR);

% Check for new indices
% Extract indices of the corresponding signal resource type
[~,sigResIndicesInfo] = extractPHYSignalResourceTypeIndices(carrier,phySigRes);
L2SMIndicesInfo = L2SMConfig.IndicesInfo;
isNewIndices = isempty(L2SMIndicesInfo) || ~isequal(L2SMIndicesInfo,sigResIndicesInfo);

% Check for new xOverhead
L2SMXOverhead = L2SMConfig.CQI.XOverhead;
isNewXOverhead = isempty(L2SMConfig.CQI.XOverhead) || ~isequal(L2SMXOverhead,xOverhead);

% Check for overall change
hasSomethingChanged = isNewTable || isNewIndices || isNewXOverhead;

if hasSomethingChanged
    % Update the L2SM configuration cache
    L2SMConfig.CQI.Table = selectedCQIQmTCR;
    L2SMConfig.CQI.XOverhead = xOverhead;
     
    % G is the physical channel bit capacities cache
    % Note: G is not Shannon capacity but is mutual information/efficieny
    % lookup values that the L2S mapping interface uses under the hood
    % In other words, G is precomputed bit "capacities" that apply to any
    % physical channels
    L2SMConfig.CQI.G = [];

    % Physical channel indices information
    L2SMConfig.IndicesInfo = sigResIndicesInfo;

    % Combine all CQI indices w.r.t codewords
    L2SMConfig = combineCQICW(L2SMConfig);
end

% Compute the effective SINR, code rate, code block  BLER and the number of 
% code block. Then map those parameters to each CQI combination.
[L2SMConfig,effSINR,codeBLER] = L2SMapping.ComputeCQISINRCodeBLER(L2SMConfig,...
    phySigRes,sinr);

% Calculate transport BLER
numCodeBlock = L2SMConfig.CQI.C;

% Probability that a single CB is correct
P_CB_correct = 1 - codeBLER;

% Probability that CBs in the TB are correct (assuming independent CB
% errors)
P_allCB_correct = P_CB_correct.^numCodeBlock;

% TB has at least one CB error -> TB fails
transportBLER = 1 - P_allCB_correct;

% Select the CQI combination with the largest BLER less than or
% equal to the threshold 'blerThresh'
idxThresh = find(all(transportBLER <= blerThresh,2),1,'last');

% If no CQI combination meets BLER criterion, select the first CQI
% combination and increase the CQI for each codeword so as to meet the BLER
% criterion
if isempty(idxThresh)
    idxThresh = checkIncreaseCQIBLERThreshold(L2SMConfig,transportBLER,blerThresh);
end
tableRow = L2SMConfig.CQI.TableRowCombos(idxThresh,:);
cqiIdx = tableRow - 1;
% Create CQI info structure
cqiInfo.EffectiveSINR = effSINR(idxThresh,:);
cqiInfo.TransportBlockSize = L2SMConfig.CQI.TransportBlockSize(idxThresh,:);
cqiInfo.Qm = L2SMConfig.CQI.Table(tableRow,1);
cqiInfo.TargetCodeRate = L2SMConfig.CQI.Table(tableRow,2) / 1024;
cqiInfo.G = L2SMConfig.CQI.G(idxThresh,:);
cqiInfo.NBuffer = L2SMConfig.CQI.NBuffer(idxThresh,:);
cqiInfo.EffectiveCodeRate = L2SMConfig.CQI.EffectiveCodeRate(idxThresh,:);
cqiInfo.CodeBLER = codeBLER(idxThresh,:);
cqiInfo.C = numCodeBlock(idxThresh,:);
cqiInfo.TransportBLER = transportBLER(idxThresh,:);
end

%% Local Helper Fcn
function [ind,indInfo] = extractPHYSignalResourceTypeIndices(carrier,sigResource)

sigResType = class(sigResource);
switch sigResType
    case 'nrPDSCHConfig' % If is PDSCH
        [ind,indInfo] = nrPDSCHIndices(carrier,sigResource);
    case 'nrPUSCHConfig' % If is PUSCH
        [ind,indInfo] = nrPUSCHIndices(carrier,sigResource);
    case 'nrCSIRSConfig' % If is CSI-RS
        [ind,indInfo] = nrCSIRSIndices(carrier,sigResource);
    case 'nrSRSConfig'   % If is SRS
        [ind,indInfo] = nrSRSIndices(carrier,sigResource);
    otherwise
        error('Unknown physical signal resource type')
end
end

function L2SMConfig = combineCQICW(L2SMConfig)

cqiTable = L2SMConfig.CQI.Table(:,1);

% Find how many distinct Qm values from the selected CQI Table
QmVal = unique(cqiTable);
QmVal(isnan(QmVal)) = [];

% Available modulation scheme (constellation mapping) and its mod order
validModScheme = {'pi/2-BPSK'; 'QPSK'; '16QAM'; '64QAM'; '256QAM'; '1024QAM'};
validModOrd    = [          1;      2;       4;       6;        8;       10]; 

% Find its respective mod scheme to QmVal
selectedModScheme = cell(length(QmVal),1);
for b = 1:length(QmVal)
    selectedIdx = validModOrd == QmVal(b);
    selectedModScheme{b} = validModScheme{selectedIdx};
end

% Get CQI Table row index combination
numCQI = size(cqiTable);
numCodewords = numel(L2SMConfig.IndicesInfo.G);
tableTotalSize = numCQI(1)^numCodewords;
tableRowCombos = ind2sub(numCQI,(1:tableTotalSize).');

% Update CQI combinations and QmVal
L2SMConfig.CQI.TableRowCombos = tableRowCombos;
L2SMConfig.CQI.TableQmValues = QmVal;
L2SMConfig.CQI.TableModulations = selectedModScheme;
end

function bestIdx = checkIncreaseCQIBLERThreshold(L2SMConfig,transportBLER,blerThresh)

% First index always selected
bestIdx = 1;

numCodewords = size(transportBLER,2);
if numCodewords > 1

    % Per codeword
    cwIdxs = 1:numCodewords;
    for currcwIdx = cwIdxs

        % Get the current CQI combo
        thisCQIRow = L2SMConfig.CQI.TableRowCombos(bestIdx,:);

        % In multi-codeword (MCW) transmission, you can have 1 or 2 
        % transport blocks (codewords) per UE.
        % Each codeword can, in principle, use a different CQI (since 
        % it may see different SINR, precoding, etc.).
        %
        % So evaluate all possible CQI assignments across all codewords. 
        % But when evaluating one particular codeword, "hold constant" 
        % the CQI choices of the other codewords and vary only the current one.

        % For the current codeword, collect all CQI table rows where this CW's CQI
        % varies across all options, while the other CWs keep their current CQIs.

        % Mask of all codeword indices except the current one
        maskOtherCWs = (cwIdxs ~= currcwIdx);

        % Extract the CQIs of the other codewords from all rows of the combo table
        otherCWsAllRows = L2SMConfig.CQI.TableRowCombos(:,maskOtherCWs);

        % Extract the CQIs of the other codewords from the current table row
        otherCWsCurrentRow = thisCQIRow(maskOtherCWs);

        % Compare all rows' "other CWs CQIs" against the current row's values
        matchMatrix = otherCWsAllRows == otherCWsCurrentRow;

        % Keep only rows where all the other CWs match (row-wise AND)
        matchingRows = all(matchMatrix,2);

        % Get indices of those rows in the combo table
        allIdxCQIsThisCW = find(matchingRows);

        % Then find the row that meets the BLER criterion for the current
        % codeword
        criterionMet = transportBLER(allIdxCQIsThisCW,currcwIdx) <= blerThresh;
        bestIdxCQIThisCW = find(criterionMet,1,'last');

        if ~isempty(bestIdxCQIThisCW)
            bestIdx = allIdxCQIsThisCW(bestIdxCQIThisCW);
        end
    end
end
end