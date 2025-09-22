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
L2SMXOverhead = L2SM.CQI.XOverhead;
isNewXOverhead = isempty(L2SM.CQI.XOverhead) || ~isequal(L2SMXOverhead,xOverhead);

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
[L2SMConfig,effSINR,codeBLER] = ComputeCQISINRCodeBLER(L2SMConfig,phySigRes,sinr);

% Calculate transport BLER
numCodeBlock = L2SMConfig.CQI.C;

% Probability that a single CB is correct
P_CB_correct = 1 - codeBLER;

% Probability that CBs in the TB are correct (assuming independent CB
% errors)
P_allCB_correct = P_CB_correct^numCodeBlock;

% TB has at least one CB error -> TB fails
transportBLER = 1 - P_allCB_correct;

% Select the CQI combination with the largest BLER less than or
% equal to the threshold 'blerThresh'
idxThresh = find(all(transportBLER <= blerThresh,2),1,'last');

end

function [ind,indInfo] = extractPHYSignalResourceTypeIndices(carrier,sigResource)

sigResType = class(sigResource);
switch class(sigResType)
    case 'nrPDSCHConfig' % If is PDSCH
        [ind,indInfo] = nrPDSCHIndices(carrier,sigResType);
    case 'nrPUSCHConfig' % If is PUSCH
        [ind,indInfo] = nrPUSCHIndices(carrier,sigResType);
    case 'nrCSIRSConfig' % If is CSI-RS
        [ind,indInfo] = nrCSIRSIndices(carrier,sigResType);
    case 'nrSRSConfig'   % If is SRS
        [ind,indInfo] = nrSRSIndices(carrier,sigResType);
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
validModScheme = {'pi/2-BPSK', 'QPSK', '16QAM', '64QAM', '256QAM', '1024QAM'};
validModOrd    = [          1,      2,       4,       6,        8,       10]; 

% Find its respective mod scheme to QmVal
selectedModScheme = cell(length(QmVal),1);
for b = 1:length(QmVal)
    selectedIdx = validModOrd == QmVal;
    selectedModScheme{b} = validModScheme{selectedIdx};
end

% Get CQI Table row index combination
numCQI = size(cqiTable);
numCodewords = numel(L2SMConfig.IndicesInfo.G);
tableTotalSize = numCQI^numCodewords;
[r,c] = ind2sub(numCQI,(1:tableTotalSize).');
tableRowCombos = [r,c];

% Update CQI combinations and QmVal
L2SMConfig.CQI.TableRowCombos = tableRowCombos;
L2SMConfig.CQI.TableQmValues = QmVal;
L2SMConfig.CQI.TableModulations = selectedModScheme;
end

function bestPerCW = pickBestPerCodeword()

end