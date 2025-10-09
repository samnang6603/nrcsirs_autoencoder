function [modulation,targetCodeRate,precoder] = DecodeCSI(carrier,pdsch,...
    pdschX,csiReport,csiFeedbackOpts)
%DecodeCSI decodes CSI to obtain new modulation, target code rate and
%precoding matrix


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

subbandInfo = CSIReporting.GetDLPMISubbandInfo(reportConfig);
wtx = CSIReporting.ComputePRGPrecoders(carrier,prgBundleSize,W,subbandInfo.SubbandSet);

end
