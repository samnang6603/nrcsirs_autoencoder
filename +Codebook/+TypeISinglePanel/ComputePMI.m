function codebook = ComputePMI(reportConfig,nHEstCSIRSPort,numLayers)
%COMPUTEPMI Calculate PMI and get the codebook for Type I single-panel
%codebook.

codebookMode = reportConfig.CodebookMode;
codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
i2Restriction = reportConfig.i2Restriction;

% Co-phasing factor value fcn handle
% TS 38.214 5.2.2.2.1
phi = @(x)exp(1i*pi*x/2);

if nHEstCSIRSPort == 2 % basically simple 2 Tx
    % For numLayer = 1 & 2, CSI reporting uses ports 3000 - 3001 according
    % to TS 38.214 Table 5.2.2.2.1-1
    if numLayers == 1
        % Means take the two Tx and combine them with different relative
        % phases, yielding beam steering in different directions
        codebook = zeros(0,0,4);
        codebook(:,:,1) = 1/sqrt(2).*[1;  1 ]; % no phase offset
        codebook(:,:,2) = 1/sqrt(2).*[1;  1i]; % 90 deg phase offset
        codebook(:,:,3) = 1/sqrt(2).*[1; -1 ]; % 180 deg phase offset
        codebook(:,:,4) = 1/sqrt(2).*[1; -1i]; % 270 deg phase offset
        
        % Grab indices of precoders forbidden by subset restriction
        restrictedIndices = find(~codebookSubsetRestriction);

        % For numLayers = 1, only the first 4 entries of the codebook are
        % valid (defined above as 1-4 in linear indices)
        restrictedIndices = restrictedIndices(restrictedIndices <= 4);
        if ~isempty(restrictedIndices)
            % if no restriction indices, check if 1-4 were restricted
            % 1,2,3,4 are the linear indices of the codebook entries
            %restrictedSet = sum(restrictedIndices == [1;2;3;4],2);
            restrictedSet = intersect(restrictedIndices,1:4);

            % if there were any restricted, zero them out
            codebook(:,:,restrictedSet) = 0;
        end
    elseif numLayers == 2
        % If numLayers = 2, then use the 2nd column of Table 5.2.2.2.1-1
        codebook(:,:,1) = 1/2*[1 1;1 -1]; % orthogonal matrix no offset
        codebook(:,:,2) = 1/2*[1 1; 1i -1i]; % orthogonal matrix 90 deg offset

        % Same logic as numLayer = 1
        % For numLayers = 2, only the 5th and 6th entries of the codebook
        % are valid (defined above as 1-2 in linear indices)
        restrictedIndices = find(~codebookSubsetRestriction);
        restrictedIndices = restrictedIndices(restrictedIndices > 4);
        if ~isempty(restrictedIndices)
            % if no restriction indices, check if 5 & 6 were restricted
            % 5,6 are the linear indices of the codebook entries,
            % continuing from 1,2,3,4 for numLayer = 1 case
            %restrictedSet = sum(restrictedIndices == [5;6],2);
            restrictedSet = intersect(restrictedIndices,[5,6]);
            codebook(:,:,restrictedSet) = 0;
        end
    end

elseif nHEstCSIRSPort > 2
    % For more than 2 Tx
    panelDim = reportConfig.PanelDimensions;
    N1 = panelDim(1);
    N2 = panelDim(2);
    O1 = reportConfig.OverSamplingFactors(1);
    O2 = reportConfig.OverSamplingFactors(2);

    if numLayers == 1
        % Codebooks for numLayers = 1 uses antenna ports 
        % 3000 - 2999+nHEstCSIRSPort according to TS 38.214 Table
        % 5.2.2.2.1-5
        if codebookMode == 1
            i11Len = N1*O1;
            i12Len = N2*O2;
            i2Len = 4;
            codebook = zeros(nHEstCSIRSPort,numLayers,i2Len,i11Len,i12Len);
            % Exhuast all permutations of i11, i12 and i2
            for i11 = 0:i11Len-1
                for i12 = 0:i12Len-1
                    for i2 = 0:i2Len-1
                        l = i11;
                        m = i12;
                        n = i2;
                        strideIdx = N2*O2*l+m;
                    end
                end
            end

        end


    elseif numLayers == 2

    elseif (numLayers == 3) || (numLayers == 4)

    elseif (numLayers == 5) || (numLayers == 6)

    else % numLayers = 7 or 8

    end



end