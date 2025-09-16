function [PMISet,PMIInfo] = DLPMISelect(carrier,csirs,csirsIndSub,reportConfig,numLayers,H,nVar)
%DLPMISelect PMI Calculation for Downlink (PDSCH)

% Source of values: 
%   - 3GPP TS 38.214 Table 5.2.2.2.1-2 (Type I single-panel codebook)
% Inspired by:
%   - MathWorks 5G Toolbox (hDLPMISelect.m), refactored into table form 
%     for clarity, maintainability, and bug prevention.
%

reportConfig = Codebook.TypeISinglePanel.ConfigureCodebookParameters(reportConfig);

[codebook,codebookIdxSetSizes] = getCodebook(reportConfig,csirs.NumCSIRSPorts,numLayers);


end



function [codebook,codebookIdxSetSizes] = getCodebook(reportConfig,numCSIRSPorts,numLayers)
%
%

switch reportConfig.CodebookType
    case 'TypeISinglePanel'
        if numCSIRSPorts == 1
            % If the number of port is 1, then codebook is a scalar value
            % of unity
            codebook = 1;
        else
            % Codebook is an N-dim of size 
            % nHEstCSIRSPorts x numLayers x i2Len x i11Len x i12Len x i13Len
            % - OR -
            % nHEstCSIRSPorts x numLayers x i2Len
            reportConfig = Codebook.TypeISinglePanel.ConfigureCodebookParameters();
            [codebook,codebookSize] = Codebook.TypeISinglePanel.ComputePMI(reportConfig,numCSIRSPorts,numLayers);
        end

        % Size of the codebook
        i2Len  = codebookSize(3);
        i11Len = codebookSize(4);
        i12Len = codebookSize(5);
        i13Len = codebookSize(6);

        codebookIdxSetSizes = [i2Len; i11Len; i12Len; i13Len];

    case 'TypeIMultiPanel'

    case 'TypeII'

end
end
