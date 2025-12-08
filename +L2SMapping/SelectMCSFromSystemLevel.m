function [mcsIdx,mcsTbRow,transportBLER] = SelectMCSFromSystemLevel(...
    carrier,pdsch,xOverhead,SINRs,mcsTableName,blerThreshold)
%SelectMCS select MCS based on system level BLER and SINR calculations

% Intialize L2S Mapping structure
L2SMConfig = L2SMapping.Initialize(carrier);

% Number of codewords
numCodewords = pdsch.NumCodewords;

% Allocate output as NaN
mcsIdx = NaN(1,numCodewords);
mcsTbRow = NaN(1,4); % Max 4 rows
transportBLER = NaN(1,numCodewords);

% SINR per layer for real numbers
SINRs = 10*log10(SINRs + eps(SINRs));
nanIdx = isnan(SINRs);
% If there are at least one NaN across all layers, then return to invoking 
% function with the allocated outputs
if any(nanIdx,'all')
    return
end

% Extract the SINR from the NaN indices
SINRs = SINRs(~any(nanIdx),:);

% Extract modulation orders and target code rates from MCS table
mcsTable = extractMCSTable(mcsTableName);

% MCS selection from System level BLER threshold and effective SINRs
Qm = mcsTable(:,2);
TCR = mcsTable(:,3);
[L2SMConfig,mcsIdx,mcsInfo] = L2SMapping.MapSINR2CQI(L2SMConfig,carrier,...
    pdsch,xOverhead,SINRs,Qm,TCR,blerThreshold);

% Output the rest from mcsIdx
mcsTbRow = mcsTable(mcsIdx+1,:);
transportBLER = mcsInfo.TransportBLER;

end

function mcsTb = extractMCSTable(tbName)
%extractMCSTable extracts MCS table from input table name

persistent tb; % persistent variable to avoid reloading every iteration

if isempty(tb)
    mcsTbClass = nrPDSCHMCSTables;
    alltb = ["QAM64Table","QAM256Table","QAM64LowSETable","QAM1024Table"];
    for tbIdx = 1:length(alltb)
        tmpTb = mcsTbClass.(alltb(tbIdx));
        lut = [             tmpTb.MCSIndex,...
                                  tmpTb.Qm,...
               (tmpTb.TargetCodeRate)*1024,...
                   tmpTb.SpectralEfficieny];

        % Select from LUT array
        tcr = lut(:,3);
        tb{tbIdx} = lut(~isnan(tcr),:);
    end
end
tbNames = ["Table1","Table2","Table3","Table4"];
mcsTb = tb{tbName == tbNames};
end