function [RI,PMISet] = RISelect(carrier,csirs,dmrsConfig,reportConfig,H,varargin)
%RISelect Calculates and returns RI and PMISet
%   Detailed explanation goes here

% Start with a small noise var
nVar = 1e-10;
alg = 'MaxSINR'; % default algorithm
if nargin > 5
    nVar = varargin{1}; % This argument in the noise var
elseif nargin > 6
    nVar = varargin{1};
    alg = varargin{2};
end

% Get Codebook Type
cbType = reportConfig.CodebookType;

% Get CSI-RS indices
csirsInd = nrCSIRSIndices(carrier,csirs,IndexStyle="subscript");

% Update BWP and Total RBs from carrier
reportConfig.NStartBWP = carrier.NStartGrid;
reportConfig.NSizeBWP = carrier.NSizeGrid;

% Make RI restriction vector, each element corresponds to each rank (1-8)
RIRestrictionVector = makeRIVector(reportConfig);
reportConfig.RIRestriction = RIRestrictionVector;

% Calculate the number of subbands and size of each given configuration
PMISubbandInfo = CSIReport.DLPMISubbandInfo(reportConfig);

% Number of CSI-RS port
nPortCSIRS = size(H,4);

% Number of RX
nRxAnts = size(H,3);

% Find max possible rank defined in TS 38.214 5.2.2.2.1 
switch cbType
    case 'TypeISinglePanel'
        % Max rank is 8 for this type
        maxRank = min([8,nPortCSIRS,nRxAnts]); % if RX or Port exceeds 8
    case 'TypeII'
        maxRank = min(nRxAnts,2);
    case {'ETypeII','TypeIMultipanel'}
        maxRank = min(nRxAnts,4);
    otherwise
        maxRank = min(nRxAnts,4);
end

% Find unrestricted ranks that are less than or equal to max rank
unrestrictedRanks = find(reportConfig.RIRestriction);
possibleRanks = intersect(1:maxRank,unrestrictedRanks);

% Initialize RI and PMISet allocation
RI = [];
PMISet = struct();
switch cbType
    case 'Type1SinglePanel'
        PMISet.i1 = [NaN, NaN, NaN];
        PMISet.i2 = Nan(1,PMISubbandInfo.NumSubbands);
    otherwise
        PMISet.i1 = NaN(1,6);
        PMISet.i2 = Nan(3,PMISubbandInfo,NumSubbands);
end

CQIPMICalcParams = struct();
CQIPMICalcParams.Carrier = carrier;
CQIPMICalcParams.CSIRS = csirs;
CQIPMICalcParams.DMRSConfig = dmrsConfig;
CQIPMICalcParams.ReportConfig = reportConfig;
CQIPMICalcParams.Channel = H;
CQIPMICalcParams.NoiseVariance = nVar;
CQIPMICalcParams.PossibleRanks = possibleRanks;
CQIPMICalcParams.PMISubbandInfo = PMISubbandInfo;

if ~isempty(possibleRanks) && ~isempty(csirsInd)
    if strcmpi(alg,'MaxSE')
        [RI,PMISet] = calculateCQI(CQIPMICalcParams,PMISet);
    else
        [RI,PMISet] = calculatePMI(CQIPMICalcParams,PMISet);
    end
end



end

function [RI,PMISet] = calculateCQI(CQIPMICalcParams)
%calculateCQI - CQI Computation

possibleRanks = CQIPMICalcParams.PossibleRanks;


% Extract spectral efficiency from standard CQI table TS 38.211
cqiTableName = ['CQI' CQIPMICalcParams.ReportConfig.CQITable];
cqiTable = nrCQITables;
spectralEfficiency = cqiTable.(cqiTableName).SpectralEfficiency;

% Find the best CQI for each possible rank, then select the rank that
% yields highest coding and modulation efficiencies
maxRank = max(possibleRanks);
pmi = 1:maxRank;
efficiency = NaN(maxRank,1);
for r = possibleRanks
    % Find CQI and PMI for current rank
    CQIPMICalcParams.ThisRank = r;
    [cqi,pmi(r),cqiInfo] = s;



end

end









function RIRestrictionVector = makeRIVector(reportConfig)
% Validate and make RI vector with each element corresponds to the index of
% the individual rank, with max dictates by maxRank
switch reportConfig.CodebookType
    case 'Type1SinglePanel'
        maxRank = 8;
        %codebookType = 'type I single-panel';
    case 'Type1MultiPanel'
        maxRank = 4;
        %codebookType = 'type I multi-panel';
    case 'Type2'
        maxRank = 2;
        %codebookType = 'type II';
    case 'EType2'
        maxRank = 4;
        %codebookType = 'enhanced type II';
    otherwise
        error('Invalid Codebook Type');
end
RIRestrictionVector = ones(1,maxRank);
end