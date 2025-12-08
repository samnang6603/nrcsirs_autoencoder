function [codebook,codebookIdxSetSizes] = ExtractCodebook(reportConfig,Pcsirs,numLayers)
%COMPUTEPMI Calculate PMI and get the codebook for Type I single-panel
%codebook.

codebookMode = reportConfig.CodebookMode;
codebookSubsetRestriction = reportConfig.CodebookSubsetRestriction;
i2Restriction = reportConfig.i2Restriction;

% Co-phasing factor value fcn handle
% TS 38.214 5.2.2.2.1
%phi = @(x)exp(1i*pi*x/2);

if Pcsirs == 2 % basically simple 2 Tx
    % For numLayer = 1 & 2, CSI reporting uses ports 3000 - 3001 according
    % to TS 38.214 Table 5.2.2.2.1-1
    if numLayers == 1
        % Means take the two Tx and combine them with different relative
        % phases, yielding beam steering in different directions
        codebook = zeros(2,1,4); % Pre-allocation
        codebook(:,:,1) = 1/sqrt(2).*[1;  
                                      1]; % no phase offset
        codebook(:,:,2) = 1/sqrt(2).*[1;  
                                      1i]; % 90 deg phase offset
        codebook(:,:,3) = 1/sqrt(2).*[1; 
                                     -1]; % 180 deg phase offset
        codebook(:,:,4) = 1/sqrt(2).*[1; 
                                     -1i]; % 270 deg phase offset
        
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

elseif Pcsirs > 2
    % For more than 2 Tx
    panelDim = reportConfig.PanelDimensions;
    N1 = panelDim(1);
    N2 = panelDim(2);
    O1 = reportConfig.OverSamplingFactors(1);
    O2 = reportConfig.OverSamplingFactors(2);

    if numLayers == 1
        % Codebooks for numLayers = 1 uses antenna ports 
        % 3000 - 2999+Pcsirs according to TS 38.214 Table
        % 5.2.2.2.1-5
        if codebookMode == 1
            i11Len = N1*O1;
            i12Len = N2*O2;
            i2Len = 4;
            codebook = zeros(Pcsirs,numLayers,i2Len,i11Len,i12Len);
            % Exhuast all permutations of i11, i12 and i2
            for i11 = 0:i11Len-1
                for i12 = 0:i12Len-1
                    for i2 = 0:i2Len-1
                        l = i11;
                        m = i12;
                        n = i2;
                        % linStrideIdx:
                        % Flattened index = i11 * (N2*O2) + i12
                        % Flattened index = i11*i12Len + i12
                        % Each row (i11) is a stride of N2*O2, i12 offsets within the stride.
                        linStrideIdx = l*N2*O2*+m;
                        [lmRestricted,i2Restricted] = checkRestriction(codebookSubsetRestriction,linStrideIdx,n,i2Restriction);
                        if ~(lmRestricted || i2Restricted)
                            vlm = computeVlm(N1,N2,O1,O2,l,m);
                            phi_n = exp(1i*pi*n/2); % co-phasing factor
                            % Calculate Wlmn from Table 5.2.2.2.1-5
                            tmp = [       vlm; 
                                    phi_n*vlm];
                            W1 = 1/sqrt(Pcsirs)*tmp;
                            codebook(:,:,n+1,l+1,m+1) = W1;
                        end
                    end
                end
            end
        else % codebookMode == 2
            % To be done later

        end


    elseif numLayers == 2
        % numLayers == 2 uses TS 38.214 Table 5.2.2.2.1-6 via antenna ports
        % 3000 - 2999+Pcsirs

        % numLayers = 2 requires i13, k1 and k2
        % Using TS 38.214 Table 5.2.2.2.1-3
        % The values are arrange via index i13 = [0,1,2,3]
        if (N1 > N2) && (N2 > 1)
            k1 = [0, O1,  0, 2*O1];
            k2 = [0,  0, O2,    0];
            i13Len = 4;
        elseif (N1 == N2)
            k1 = [0, O1,  0,   O1];
            k2 = [0,  0, O2,   O2];
            i13Len = 4;
        elseif (N1 == 2) && (N2 == 1)
            k1 = [0, O1];
            k2 = [0,  0];
            i13Len = 2;
        else
            k1 = [0, O1, 2*O1, 3*O2];
            k2 = [0,  0,    0,    0];
            i13Len = 4;
        end

        if codebookMode == 1
            i11Len = N1*O1;
            i12Len = N2*O2;
            i2Len = 2;
            codebook = zeros(Pcsirs,numLayers,i2Len,i11Len,i12Len,i13Len);
            % Exhuast all permutations of i11, i12, i13 and i2
            for i11 = 0:i11Len-1
                for i12 = 0:i12Len-1
                    for i13 = 0:i13Len-1
                        for i2 = 0:i2Len-1
                            l  = i11;
                            lp = i11 + k1(i13+1); % k1 is the mini LUT, i3 is the stride
                            m  = i12;
                            mp = i12 + k2(i13+1);
                            o  = i13;
                            n  = i2;
                            linStrideIdx = l*N2*O2 + m;
                            [lmRestricted,i2Restricted] = checkRestriction(codebookSubsetRestriction,linStrideIdx,n,i2Restriction);
                            if ~(lmRestricted || i2Restricted)
                                vlm = computeVlm(N1,N2,O1,O2,l,m);
                                vlpmp = computeVlm(N1,N2,O1,O2,lp,mp); % v l_prime p_prime
                                phi_n = exp(1i*pi*n/2);
                                tmp = [      vlm,        vlpmp; 
                                       phi_n*vlm, -phi_n*vlpmp];
                                W2 = (1/sqrt(2*Pcsirs))*tmp;
                                codebook(:,:,n+1,l+1,m+1,o+1) = W2;
                            end
                        end
                    end
                end
            end
        end

    elseif (numLayers == 3) || (numLayers == 4)
        % Combine calculation because the only difference between 3 and 4 
        % layers is the dimension of the precoder W. The procedures are
        % identical.

        if Pcsirs < 16
            % If PCSIRS < 16 then use TS 38.214 Table 5.2.2.2.1-4 to
            % compute i13 indices which corresponds to k1 and k2.
            if (N1 == 2) && (N2 == 1)
                i13Len = 1;
                k1 = O1;
                k2 = 0;
            elseif (N1 == 4) && (N2 == 1)
                i13Len = 3;
                k1 = O1*(1:3);
                k2 = [0, 0, 0];
            elseif (N1 == 6) && (N2 == 1)
                i13Len = 4;
                k1 = O1*(1:4);
                k2 = [0, 0, 0, 0];
            elseif (N1 == 2) && (N2 == 2)
                i13Len = 3;
                k1 = [O1,  0, O1];
                k2 = [ 0, O2, O2];
            elseif (N1 == 3) && (N2 == 2)
                i13Len = 4;
                k1 = [O1,  0, O1, 2*O1];
                k2 = [ 0, O2, O2,    0];
            end
            i11Len = N1*O1;
            i12Len = N2*O2;
            i2Len  = 2;
            codebook = zeros(Pcsirs,numLayers,i2Len,i11Len,i12Len,i13Len);
            % Exhuast all permutations of i11, i12, i13, and i2
            for i11 = 0:i11Len-1
                for i12 = 0:i12Len-1
                    for i13 = 0:i13Len-1
                        for i2 = 0:i2Len-1
                            l  = i11;
                            lp = i11 + k1(i13+1);
                            m  = i12;
                            mp = i12 + k2(i13+1);
                            n  = i2;
                            o  = i13;
                            linStrideIdx = m + N2*O2*l;
                            [lmRestricted,i2Restricted] = checkRestriction(codebookSubsetRestriction,linStrideIdx,n,i2Restriction);
                            if ~(lmRestricted || i2Restricted)
                                vlm = computeVlm(N1,N2,O1,O2,l,m);
                                vlpmp = computeVlm(N1,N2,O1,O2,lp,mp); % v l_prime p_prime
                                phi_n = exp(1i*pi*n/2);
                                phi_nvlm = phi_n*vlm;
                                phi_nvlpmp = phi_n*vlpmp;
                                if numLayers == 3
                                    % 3-layers CSI reporting via antenna
                                    % ports 3000 to 2999+Pcsirs in
                                    % TS 38.214 Table 5.2.2.2.1-7
                                    tmp = [     vlm,      vlpmp,       vlm;
                                           phi_nvlm, phi_nvlpmp, -phi_nvlm];
                                    W3  = (1/sqrt(3*Pcsirs))*tmp;
                                    codebook(:,:,n+1,l+1,m+1,o+1) = W3;
                                else
                                    % 4-layers CSI reporting via antenna
                                    % ports 3000 to 2999+Pcsirs
                                    % TS 38.214 Table 5.2.2.2.1-8
                                    tmp = [     vlm,      vlpmp,       vlm,       vlpmp;
                                           phi_nvlm, phi_nvlpmp, -phi_nvlm, -phi_nvlpmp];
                                    W4  = (1/sqrt(4*Pcsirs))*tmp;
                                    codebook(:,:,n+1,l+1,m+1,o+1) = W4;
                                end
                            end
                        end
                    end
                end
            end
        else
            % Pcsirs > 16
        end

    %elseif (numLayers == 5) || (numLayers == 6)

    %else % numLayers = 7 or 8

    end

end
% Size of the codebook. Ignore the first 2 dims, and extend as far
% as i13 index. Collate it in codebookIdxSetSizes
[~,~,i2Len,i11Len,i12Len,i13Len] = size(codebook);
codebookIdxSetSizes = [i2Len, i11Len, i12Len, i13Len];
end






function [vlmRestricted,i2Restricted] = checkRestriction(codebookSubsetRestriction,...
    linearStrideIdx,coPhasingFactorIdx,i2Restriction)

% Find the restricted set and set to -1
restrictedIdx = find(~codebookSubsetRestriction)-1;

% vlm restriction set to 0 (false)
vlmRestricted = 0;

% If linear stride index has one of the restricted index then vlm
% restriction is true
if ~isempty(intersect(restrictedIdx,linearStrideIdx))
    vlmRestricted = true;
end

% Now find i2 restricted list
i2RestrictedIdx = find(~i2Restriction)-1;
i2Restricted = 0;
% If the precoding matrix based on vlm and vbarlm are restricted, update
% the i2 restriction status
if any(restrictedIdx == coPhasingFactorIdx)
    i2Restricted = true;
end
end

function vlm = computeVlm(N1,N2,O1,O2,l,m)
% TS 38.214 page 47
ul = exp(1i*2*pi*l*(0:N1-1)/(O1*N1)); % horizontal
um = exp(1i*2*pi*m*(0:N2-1)/(O2*N2)); % vertical
vlm = ul.'*um; % Outer product
vlm = vlm.';
vlm = vlm(:);
end