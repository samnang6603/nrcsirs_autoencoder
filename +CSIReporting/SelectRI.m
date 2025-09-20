function [RI,PMISet] = SelectRI(carrier,csirs,csiFeedbackOpts,H,varargin)
%RISelect Calculates and returns RI and PMISet
%   Detailed explanation goes here

% Start with a small noise var
nVar = 1e-10;
alg = 'MaxSINR'; % default algorithm
if nargin == 5
    nVar = varargin{1}; % This argument in the noise var
elseif nargin > 5
    nVar = varargin{1};
    alg = varargin{2};
end

% CSI reportConfig
reportConfig = csiFeedbackOpts.CSIReportConfig;

% DMRS
dmrsConfig = csiFeedbackOpts.DMRSConfig;

% Get CSI-RS indices
csirsInd = nrCSIRSIndices(carrier,csirs,IndexStyle='subscript');

% Update BWP and Total RBs from carrier
reportConfig.NStartBWP = carrier.NStartGrid;
reportConfig.NSizeBWP = carrier.NSizeGrid;

% Make RI restriction vector, each element corresponds to each rank (1-8)
RIRestrictionVector = makeRIVector(reportConfig);
reportConfig.RIRestriction = RIRestrictionVector;

% Calculate the number of subbands and size of each given configuration
PMISubbandInfo = CSIReporting.GetDLPMISubbandInfo(reportConfig);

% Number of CSI-RS port
Pcsirs = size(H,4);

% Number of RX
nRxAnts = size(H,3);

% Find max possible rank defined in TS 38.214 5.2.2.2.1 
switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        % Max rank is 8 for this type
        maxRank = min([8,Pcsirs,nRxAnts]); % if RX or Port exceeds 8
    case 'TypeII'
        maxRank = min(nRxAnts,2);
    case {'ETypeII','eTypeII','EnhTypeII','TypeIMultipanel'}
        maxRank = min(nRxAnts,4);
    otherwise
        maxRank = min(nRxAnts,4);
end

% Find unrestricted ranks that are less than or equal to max rank
unrestrictedRanks = find(reportConfig.RIRestriction);
possibleRanks = intersect(1:maxRank,unrestrictedRanks);

% Initialize RI and PMISet allocation
RI = NaN;
PMISet = struct();
switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        PMISet.i1 = [NaN, NaN, NaN];
        PMISet.i2 = NaN(1,PMISubbandInfo.NumSubbands);
    otherwise
        PMISet.i1 = NaN(1,6);
        PMISet.i2 = NaN(3,PMISubbandInfo,NumSubbands);
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
    switch alg
        case 'MaxSE'
            [RI,PMISet] = riCalculateCQI(CQIPMICalcParams,PMISet);
        case 'MaxSINR'
            [RI,PMISet] = riCalculatePMI(CQIPMICalcParams,PMISet);
        otherwise
            error('Unknown Criteria: Must be MaxSE or MaxSINR');
    end
end
end

function [RI,PMISet] = riCalculateCQI(CQIPMICalcParams,PMISet)
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
    [cqi,pmi(r),cqiInfo] = CSIReporting.SelectCQI(CQIPMICalcParams);



end

end




function RIRestrictionVector = makeRIVector(reportConfig)
% Validate and make RI vector with each element corresponds to the index of
% the individual rank, with max dictates by maxRank
switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        maxRank = 8;
        %codebookType = 'type I single-panel';
    case 'TypeIMultiPanel'
        maxRank = 4;
        %codebookType = 'type I multi-panel';
    case 'TypeII'
        maxRank = 2;
        %codebookType = 'type II';
    case {'eTypeII','ETypeII','EnhTypeII'}
        maxRank = 4;
        %codebookType = 'enhanced type II';
    otherwise
        error('Invalid Codebook Type');
end
RIRestrictionVector = ones(1,maxRank);
end