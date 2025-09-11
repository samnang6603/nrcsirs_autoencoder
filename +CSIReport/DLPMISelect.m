function [PMISet,PMIInfo] = DLPMISelect(carrier,csirs,csirsIndSub,reportConfig,numLayers,H,nVar)
%DLPMISelect PMI Calculation for Downlink (PDSCH)

% Source of values: 
%   - 3GPP TS 38.214 Table 5.2.2.2.1-2 (Type I single-panel codebook)
% Inspired by:
%   - MathWorks 5G Toolbox (hDLPMISelect.m), refactored into table form 
%     for clarity, maintainability, and bug prevention.
%


end



function [codebook,indexSetSizes] = getCodebook(reportConfig,numCSIRSPorts,numLayers)
%
%

switch reportConfig.CodebookType
    case 'Type1SinglePanel'
        if numCSIRSPorts == 1
            % If the number of port is 1, then codebook is a scalar value
            % of unity
            codebook = 1;
        else
            % Codebook is an N-dim of size 
            % nHEstCSIRSPorts x numLayers x i2Len x i11Len x i12Len x i13Len
            % - OR -
            % nHEstCSIRSPorts x numLayers x i2Len
            codebook = getPMIType1SinglePanelCodebook();
        end
end
end
