function [modulation,targetCodeRate,precoder] = DecodeCSI(carrier,pdsch,...
    pdschX,csiReport,csiFeedbackOpts)
%DecodeCSI decodes CSI to obtain new modulation, target code rate and
%   precoding matrix


csiReportConfig = csiFeedbackOpts.CSIReportConfig;

switch csiReportConfig.Mode
    case 'RI-PMI-CQI'
        % Map CSI to MCS
        [modulation,targetCodeRate] = mapCSI2MCS(csiReportConfig.CQITable,csiReport);

        precoder = PMISubband2PRGPrecodingMatrix(carrier,...
            pdschX.PRGBundleSize,csiReportConfig,csiReport.Precoder);

    case 'AI CSI compression'


end


end

%% Local Helper Fcn
function [modu,tcr] = mapCSI2MCS(cqiTableNameInput,csiReport)

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

modu = thisTable.Modulation(cqi+1); % add 1 due to 0-based index
% Remove 'Out of Range' row
modu = modu(~strcmpi(modu,'Out of Range'));
tcr  = thisTable.TargetCodeRate(cqi+1);
end

function wtx = PMISubband2PRGPrecodingMatrix(carrier,prgBundleSize,reportConfig,W)
% Map codebook-based precoding matrices from subbands to PRGs

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

function [wtx,sinr] = computeSVDPrecodingMatrix(carrier,numLayers,prbSet,Hest,nVar,prgBundleSize)
%computeSVDPrecodingMatrix calculates precoding matrices for all PRGs in 
% the carrier that overlap with the PDSCH allocation

persistent pdsch
if isempty(pdsch)
    pdsch = nrPDSCHConfig;
end

pdsch.NumLayers = numLayers;
pdsch.PRBSet = prbset;



end

function sinr = computePrecodedSINRSubroutine(H,W,nVar)
% Reference
% nrPrecodedSINR(H,nVar,W) from Mathworks
% Li, Ping, et al. "On the Distribution of SINR for the MMSE MIMO Receiver 
% and Performance Analysis." IEEE Transactions on Information Theory, 
% vol. 52, no. 1, 1 Jan. 2006, pp. 271–286, 
% https://doi.org/10.1109/tit.2005.860466. Accessed 23 Sept. 2023.

% Rerrange H matrix, so H is now nRxAnt-Pcsirs-K*L
H = permute(H,[2, 3, 1]); % Toss the BWP sc to the third dim to match W

R = pagemtimes(H,W);
[~,S,V] = pagesvd(R,'econ','vector');

SSqr = S.*S;
absVSqr = abs(V).^2;

% If H is 2D, compute SINR values using W as page
if size(H,3) == 1
    diagTerm = nVar./(nVar + SSqr + eps); % eps to prevent divide by 0
    diagTermPageT = pagetranspose(diagTerm);
    sumTerm = sum(absVSqr.*diagTermPageT,2);
    msee_i = squeeze(sumTerm); % msee of i-th stream
else % If H is n-dim
    SSqrPageT = pagetranspose(SSqr); % page transpose this
    nVarVec = nVar*ones(1,size(W,2)); % create a vector to match SSqrPageT
    diagTerm = nVar./(SSqrPageT + nVarVec + eps); % eps to prevent divide by 0
    sumTerm = sum(absVSqr.*diagTerm,2);
    msee_i = permute(sumTerm,[3, 1, 2]); % msee of i-th stream
end
sinr = 1./msee_i - 1;
end