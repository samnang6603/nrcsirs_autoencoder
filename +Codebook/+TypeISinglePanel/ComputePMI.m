function codebook = ComputePMI(reportConfig,nHEstCSIRSPort,numLayers)
%COMPUTEPMI Calculate PMI and get the codebook for Type I single-panel
%codebook.

codebookMode = reportConfig.CodebookMode;
codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
i2Restriction = reportConfig.i2Restriction;

% Co-phasing factor value fcn handle
% TS 38.214 5.2.2.2.1
phi = @(x)exp(1i*pi*x/2);

if nHEstCSIRSPort == 1
    % For numLayer = 1 & 2, CSI reporting uses ports 3000 - 3001 according
    % to TS 38.214 Table 5.2.2.2.1-1
    if numLayers == 1
        codebook = zeros(0,0,4);
        codebook(:,:,1) = 1/sqrt(2).*[1;  1 ];
        codebook(:,:,2) = 1/sqrt(2).*[1;  1i];
        codebook(:,:,3) = 1/sqrt(2).*[1; -1 ];
        codebook(:,:,4) = 1/sqrt(2).*[1; -1i];
        
    elseif numLayers == 2

    end


end