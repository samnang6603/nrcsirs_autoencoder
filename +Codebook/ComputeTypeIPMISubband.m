function [PMISet,W,SINRPerREPMI,subbandSINRs] = ComputeTypeIPMISubband(...
    codebook,PMISet,SINRPerRE,subbandInfo,codebookIdxSetSizes,numLayers,k)
%SetTypeIPMISubband Computes PMI for Type I single or multipanel for
%subband. Takes results from wideband computations as well

numSubbands = subbandInfo.NumSubbands;
numCSIRSPorts = size(codebook,1);

% Allocations some outputs
W = zeros(numCSIRSPorts,numLayers,numSubbands);
SINRPerREPMI = zeros(length(k),numLayers);
subbandSINRs = NaN(numSubbands,numLayers,codebookIdxSetSizes);

% Start of subband is index 0
subbandInit = 0;

% Lower bound of BWP
lowerBound = k >= (subbandInit*12 + 1);


% Sift through the subbands
for thisSB = 1:numSubbands

    % Subband size w.r.t subband index
    subbandSize = subbandInfo.SubbandSizes(thisSB);
    
    % Upper bound requires sb idx
    upperBound = k <= ((subbandInit + subbandSize)*12 + 1);

    % Intersection between lower and upper bound
    subbandInd = lowerBound & upperBound;

    % SINR value per subband
    % Here all SINRPerRE dim after subbandInd are selected in case of 
    % higher layer and different codebookMode
    sinrValPerSubband = SINRPerRE(subbandInd,:,:,:,:,:,:,:,:,:,:);
    
    % If CSI-RS is absent in the subband
    if all(isnan(sinrValPerSubband(:)))
        % Then set i2 as NaN
        PMISet.i2(:,thisSB) = NaN;
    else
        % Take the average of the SINRPerRE of all subband and all PMI
        % indices
        subbandSINRs(thisSB,:,:,:,:,:,:,:,:) = mean(SINRPerRE(subbandInd,:,:,:,:,:,:,:,:),1);
        
        % Add all the subband SINR values for all the layers for each i2
        % index set
        i11WB = PMISet.i1(1); % WB = WideBand
        i12WB = PMISet.i1(2);
        i13WB = PMISet.i1(3);
        switch reportConfig.CodebookType
            case 'TypeISinglePanel'
                thisSBSINRSum = sum(subbandSINRs(thisSB,:,:,i11WB,i12WB,i13WB),2);
                % Set decimal to 4 to combat fluctuation of SINR values
                thisSBSINRSum = round(thisSBSINRSum,4,'decimals');

                % Find i2 index of the max SINR for this subband
                [~,PMISet.i2(thisSB)] = max(thisSBSINRSum);
                tmpi2 = PMISet.i2(thisSB);

                % Get the Precoding matrix corresponding to that index
                W(:,:,thisSB) = codebook(:,:,tmpi2,i11WB,i12WB,i13WB);

                % Get the SINR per RE PMI
                SINRPerREPMI(subbandInd,:) = SINRPerRE(subbandInd,:,tmpi2,i11WB,i12WB,i13WB);
            case 'TypeIMultiPanel'
                % To be implemented later
        end
    end
    % Find the initial position of the next subband
    subbandInit = subbandInit + subbandInfo.SubbandSizes(thisSB);
end
end
