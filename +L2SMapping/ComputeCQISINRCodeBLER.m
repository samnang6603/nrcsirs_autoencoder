function [L2SMConfig,effSINR,cbBLER] = ComputeCQISINRCodeBLER(L2SMConfig,phySigRes,SINR)
%ComputeSINR2BLER Compute effective SINR, code rate, code block BLER and
%the number of code block, all per CQI combination.

% Get dimensions of the combination table, CQI and codeword
tableRowCombos = L2SMConfig.CQI.TableRowCombos;
numCombosCQI   = size(tableRowCombos,1);
numCodewords   = numel(L2SMConfig.IndicesInfo.G);

% Determine if there is a change in configuration such as CQI Table, 
% physical channel indices and xOverhead.
isConfigChanged = isempty(L2SMConfig.CQI.G);

% Allocate temp NaN array for convenience
tmpNaN = NaN(numCombosCQI,numCodewords);

% If there is a change, reinitialize cache entries
if isConfigChanged
    % Layers per codeword NL
    numLayersPerCodeword = floor((phySigRes.NumLayers + (0:numCodewords-1))/numCodewords);
    
    numPRB = numel(phySigRes.PRBSet);
    L2SMConfig.CQI.Qm = tmpNaN; % Qm
    L2SMConfig.CQI.G = tmpNaN;  % bit capacities
    L2SMConfig.CQI.TransportBlockSize = tmpNaN;  % TBS
    L2SMConfig.CQI.C = tmpNaN;  % number of code blocks
    L2SMConfig.CQI.NBuffer = tmpNaN; % rate matching buffer size
    L2SMConfig.CQI.EffectiveCodeRate = tmpNaN;
    dlschInfo = nrDLSCHInfo(1,0.5); % intialize dlsch info with trblk of 1 and CR of 0.5
    % Repeat dlschInfo for all cqi and codeword
    dlschInfoPerCQIPerCW = repmat(dlschInfo,numCombosCQI,numCodewords); 
    L2SMConfig.CQI.DLSCHInfo = dlschInfoPerCQIPerCW;
end

% Pre-allocate effective SINR and code block BLER
effSINR = tmpNaN;
cbBLER = tmpNaN;

% Number of RE per RB
numREPerRB = L2SMConfig.IndicesInfo.NREPerPRB;
xOverhead  = L2SMConfig.CQI.XOverhead;

% If split code block is true then make cell for effSINR and codeBlockBLER
if (L2SMConfig.SplitCodeBlock)
    effSINR = num2cell(effSINR);
    cbBLER  = num2cell(cbBLER);
end

% If cache entries are reinitialized, calculate and cache values that are
% not part of the incoming SINR
if isConfigChanged

    tableQmVal = L2SMConfig.CQI.TableQmValues; % Table of all Qm values
    sumNL = sum(numLayersPerCodeword); % total numLayers across all CW

    for c = 1:numCombosCQI

        % Get the Qm combination
        Qm = L2SMConfig.CQI.Table(tableRowCombos(c,:),1)';
        L2SMConfig.CQI.Qm(c,:) = Qm;

        % Get physical channel bit capacity 
        % Bit Capacity = ModOrder*NumUsableREPerLayer*NumLayerPerCW
        %            G =       Qm*                 Gd*numLayersPerCodeword(NL)
        % (Units) Bits = bits/sym*                sym*layers     
        L2SMConfig.CQI.G(c,:) = L2SMConfig.IndicesInfo.Gd*Qm.*numLayersPerCodeword;

        % For valid Qm values
        validQm = Qm > 0; 
        if any(validQm)
            % Extract transport block sizes and target code rate.
            % Output is assigned for all codewords due to nrTBS
            % automatically allocate number of codewords from
            % sum(numLayersPerCodeword) or sum(NL) even if its Qm is invalid.
            activeQm = Qm(validQm); % indices of active codewords
            activeQmLen = numel(activeQm);
            modulations = zeros(size(activeQm)); % preallocate
            for b = 1:activeQmLen
                qmVal = activeQm(b);
                % find corresponding modulation in the CQI table
                modulations(b) = L2SMConfig.CQI.TableModulations(tableQmVal == qmVal);
            end
            
            % Target code rates
            tcr = L2SMConfig.CQI.Table(tableRowCombos(c,:),2).'/1024; % per 1024 bits
            
            % Transport Block Size (TBS)
            tbs = nrTBS(modulations,sumNL,numPRB,numREPerRB,tcr(validQm),xOverhead);
            L2SMConfig.CQI.TransportBlockSize(c,:) = tbs;

            % DL-SCH information
            for b = 1:activeQmLen
                dlschInfo = nrDLSCHInfo(tbs(c,b),tcr(b));
                L2SMConfig.CQI.DLSCHInfo(c,b) = dlschInfo;
            end
        end
    end

    % Copy number of code blocks from DLSCHInfo into C in main CQI struct
    L2SMConfig.CQI.C = [L2SMConfig.CQI.DLSCHInfo.C];

    % Get rate matching buffer size
    % See TS 38.212 5.4.2.1
    numBitsPerCBLDPC = [L2SMConfig.CQI.DLSCHInfo.N]; % Ncb in TS 38.212
    numBitsPerCBLDPC = reshape(numBitsPerCBLDPC,size(L2SMConfig.CQI.DLSCHInfo));
    if ~isempty(L2SMConfig.Nref)
        % If Nref, limited buffer rate matching param, is present, take the
        % min of the Nref and numBitsPerCBLDPC
        numBitsPerCBLDPC = min(numBitsPerCBLDPC,L2SMConfig.Nref); 
    end
    for b = 1:length(numBitsPerCBLDPC)
        cbsInfo = L2SMConfig.CQI.DLSCHInfo(b);
        NBuffer = calculateSoftBufferSize(cbsInfo,numBitsPerCBLDPC);
        L2SMConfig.CQI.NBuffer(b) = NBuffer;
    end
    L2SMConfig.CQI.NBuffer = L2SMConfig.NBuffer.L2SMConfig.CQI.C;

    % Get effective code rate by considering the rate repetition in the
    % first RV for code rates lower than the mother code rate for the LDPC
    % base graph.
    % Take the minimum between bit capacities and the rate matching buffer
    den = min(L2SMConfig.CQI.G,L2SMConfig.CQI.NBuffer);
    L2SMConfig.CQI.EffectiveCodeRate = L2SMConfig.CQI.TransportBlockSize./den;
end

% Layer demapping on input SINR
layerSINR = nrLayerDemap(SINR);

% For each codeword, compute effective SINR calculation and code block BLER
for cwIdx = 1:numCodewords

    % --------------- Start of Effective SINR Computation -----------------
    % Gather all relevant inputs [Qm, G, NBuffer, C]
    effSINRCalcInputSets = [L2SMConfig.CQI.Qm(:,cwIdx),...
        L2SMConfig.CQI.G(:,cwIdx),L2SMConfig.CQI.NBuffer(:,cwIdx),...
        L2SMConfig.CQI.C(:,cwIdx)];
    
    % Find distinct parameter 'u' and the associated indices 'uIdx' in the
    % effSINRInputSet
    [u,~,uIdx] = unique(effSINRCalcInputSets,'rows');

    for nonnanIdx = find(~any(isnan(u),2))
        thisInputSet = u(nonnanIdx,:);

        % Calculate effective SINR for the current codeword and map the
        % results into the elements of the output that correspond to CQI
        % with parameters matching that input set.
        Qm = thisInputSet(1);
        G  = thisInputSet(2);
        NBuffer = thisInputSet(3);
        
        % Split SINR into code block segment (CBS) so the effective SINR
        % can be calcualted per code block
        if (L2SMConfig.SplitCodeBlocks)
            % Split Code Block >>>>>>>>>>>>>>> TBI >>>>>>>>>>>>>>>>>>>>>>>>
            %splitSINR = splitCodeBlocks(layerSINR(cwIdx),C); 
        else
            splitSINR = layerSINR(cwIdx);
        end

        % Calculate effective SINR considering rate repetition in the first
        % RV for code rates lower than the mother code rate for the LDPC
        % base graph
        a = L2SMConfig.Alpha;
        b = L2SMConfig.Beta;
        eTmp = effectiveSINRMapping(Qm,splitSINR,a,b,G,NBuffer);
        effSINR(uIdx==nonnanIdx,cwIdx) = eTmp;
    end
    % ----------------- End of Effective SINR Computation -----------------

    % --------------- Start of Code Block BLER Computation ----------------
    % Gather all relevant inputs [Qm, G, NBuffer, C]
    ecrTmp = round(L2SMConfig.CQI.EffectiveCodeRate(:,cwIdx)*1024); % rounded effective code rate
    codeBlockBLERCalcInputSets = [L2SMConfig.CQI.TransportBlockSize(:,cwIdx),...
        L2SMConfig.CQI.Qm(:,cwIdx), ecrTmp];

    % Find distinct parameter 'v' and the associated indices 'vIdx' in the
    % codeBlockBLERCalcInputSets
    [v,~,vIdx] = unique(codeBlockBLERCalcInputSets,'rows');

    for nonnanIdx = find(~any(isnan(v),2))
        thisInputSet = v(nonnanIdx,:);

        % Calculate code block BLER for the current codeword and map the
        % results into the elements of the output that correspond to CQI
        % with parameters matching that input set.
        trBlkSize = thisInputSet(1);
        Qm = thisInputSet(2);
        ecr = thisInputSet(3)/1024; % Effective Code Rate

        % If configuration has changed (new), cache the DL-SCH info for the
        % current input set to avoid computing it everytime 
        % computeSINRToCodeBlockBLER() is invoked (more efficient)
        if isConfigChanged
            L2SMConfig.CQI.DLSCHInfo(vIdx == nonnanIdx) = nrDLSCHInfo(trBlkSize,ecr);
        end

        % 'vIdx1' is the first index we use as a representative for this 
        % parameter set. Even though effectiveSINR may differ across 
        % subbands/codewords initially, all entries with vIdx == nonnanIdx 
        % collapse to the same final result (same effective SINR and DL-SCH 
        % info) after grouping. So we can just pick one (vIdx1).
        vIdx1 = find(vIdx == nonnanIdx,1);

        % Compute the code block BLER
        esize = size(effSINR(vIdx1,cwIdx));
        thisdlschInfo = L2SMConfig.CQI.DLSCHInfo(vIdx1,cwIdx);
        c = computeCodeBlockBLER(esize,trBlkSize,Qm,ecr,thisdlschInfo);



    end



    % --------------- End of Code Block BLER Computation --------------

end


end

function NBuffer = calculateSoftBufferSize(cbsInfo,Ncb)
%calculateSoftBufferSize Calculate soft buffer size. Procedure listed in 
% TS 38.212 5.4.2: 
% 1. Remove filler bits
% 2. Account for punctured bits
% 3. Compute soft buffer size for circular buffer

% Systematic bits puncturing
K = cbsInfo.K - 2*cbsInfo.Zc;  % total coded bits minus punctured systematic bits

% Fillers exclusion
Kd = K - cbsInfo.F; % F is number of NULL filler bits

% Number of filler bits inside the circular buffer
NFillterBits = max(min(K,Ncb)-Kd,0);

% Buffer size without filler bits
NBuffer = Ncb - NFillterBits;
end

function effSINR = effectiveSINRMapping(Qm,SINR,alph,beta,G,NBuffer)
    
if NBuffer < G
    % Adjust SINR to account for HARQ chase combining (at SLS level)
    Gsum = G;
    Csum = min(G,NBuffer);
    SINR = adjustSINRChaseCombining(SINR,Gsum,Csum);
end

numCodewords = numel(Qm);
maxConfigs  = size(SINR,1); % number of parameter sets (CQI configs etc.)
effSINR = zeros(maxConfigs,numCodewords); % Pre-allocate

% Load RBIR from 802.11-14/1450r0 - Box 0 Calibration Results, extended for
% 1024 QAM and 4096 QAM. Sourced from MATLAB. RBIR is empirical.
rbir = load('L2SMapping/rbir.mat');

for cfgIdx = 1:maxConfigs
    for cwIdx = 1:numCodewords
        inTmp = SINR{cfgIdx,cwIdx};
        numSym = 2^Qm{cwIdx};
        sinrTmp = computeEffectiveSINRSubroutine(inTmp,numSym,alph,beta,rbir);
        effSINR(cfgIdx,cwIdx) = sinrTmp;
    end
end

end

function effSINR = computeEffectiveSINRSubroutine(sinr,numSym,alph,beta,rbir)
%computeEffectiveSINRSubroutine Subroutine to compute effective SINR with
%tuning parameters alpha and beta and Received Bit Information Rate (RBIR)
%LookUp Table (LUT)
    
% Retrieve RBIR table based on modulation scheme (number of symbols)
thisIdx = numSym == rbir.tableModScheme;
rbirSNR = rbir.rbirTable{thisIdx}(:,1); % SNR at first column
rbirVal = rbir.rbirTable{thisIdx}(:,2); % Val at second column

% Tune with beta: pre-scaling SINR per RE before it gets mapped into Mutual
% Information (MI). Beta adjust how aggressively each RE contributes to MI
sinr = sinr/beta;

% Ensure the SINR values fall inside the LUT range
% If SINR too low/high, clip to table min/max avoid extrapolation
sinr(sinr < rbir.minSNR(thisIdx)) = rbir.minSNR(thisIdx);
sinr(sinr > rbir.maxSNR(thisIdx)) = rbir.maxSNR(thisIdx);

% Lookup and interpolate the RBIR for each RE based on its empirical SINR
% Basically: take input SINR and interpolate using values on an LUT that
% was built in SNR-only conditions
scrbir = interp1(rbirSNR,rbirVal,sinr); % RBIR per subcarrier
avgbir = mean(scrbir,'all'); % Take the average RBIR

% Invert the mapping: take the average RBIR and and map it back to an
% effective SINR
invMap = interp1(rbirVal,rbirSNR,avgbir);

% Tune with alpha: post-scaling SINR per RE after averaging and demapped.
% Alpha shifts the entire abstraction curve
effSINR = invMap*alph;
end

function outSINR = adjustSINRChaseCombining(inSINR,Gsum,Csum) 
%adjustSINRChaseCombining Adjust SINRs to account for the effect of HARQ
%chase combining (LLR from previous (re)transmissions).
%
% Gsum: Total number of coded bits transmitted so far across all HARQ
%       transmission for a codeword
% Csum: Normalization term; the number of unique coded bits that can be
%       stored. Bounded by minimum between Gsum and NBuffer

maxConfigs = size(inSINR,1); % number of parameter sets (CQI configs etc.)
numCodewords = size(inSINR,2);
outSINR = cell(maxConfigs,numCodewords); % pre-allocate as cell
for cfgIdx = 1:maxConfigs
    for cwIdx = 1:numCodewords
        chaseBoost = 10*log10(Gsum(cwIdx)/Csum(cwIdx)); % boost from chase combining
        outSINR{cfgIdx,cwIdx} = inSINR{cfgIdx,cwIdx} + chaseBoost; % add the boost
    end
end

end

function cBBLER = computeCodeBlockBLER(effSINR,trBlkSizes,Qm,ecr,dlschInfoList)

% Declare persistent variable to avoid having to reload this structure 
% everytime
persistent awgnTables; 
if isempty(awgnTables)
    % Load and retrieve AWGN table data
    Data  = load('L2SMapping/L2SM.mat');
    awgnTables = Data.awgnTable;
    for b = 1:size(awgnTables.BGN,1)
        for r = 1:size(awgnTables.data(b).R,1)
            awgnTables.data(b).data(r).data = double(awgnTables.data(b).data(r).data);
        end
    end
end

awgnLUT = awgnTables;

% Clamp effective code rate to be between 1/1024 and 1023/1024
ecr = max(ecr,1023/1024);
ecr = min(ecr,   1/1024);
R = round(ecr*1024); % integer code rate: rate/1024

numCodewords  = numel(trBlkSizes);
maxC = size(effSINR,1);
codeBlockBLER = zeros(effSINR,numCodewords);
% For each codeword
for cwIdx = 1:numCodewords

    % Lookup DL-SCH info
    dlschInfo = dlschInfoList(cwIdx);

    % Get BGN and Zc
    BGN = dlschInfo.BGN;
    Zc  = dlschInfo.Zc;

    % Lookup AWGN table corresponding to BGN, R, Qm, Zc
    thisLUT = awgnLUT.data(awgn.BGN == BGN); % select BGN as main index
    R_range = (R >= thisLUT.R(:,1)) & (R <= thisLUT.R(:,2));
    thisLUT = thisLUT.data(R_range);

    % Find table with the target Qm and Zc and pick the smallest available 
    % code rate in the LUT that is >= the target R (Round up conservatively
    % so the selected BLER curve does not underestimate errors.)
    iR   = find(thisLUT.R >= R,1); % find index of the diff that is 0 and greater than
    iQm  = find(thisLUT.Qm == Qm); % find index corresponding to target Qm
    iZc  = find(thisLUT.Zc == Zc); % find index corresponding to target Zc
    iLUT = thisLUT.data(:,:,iR,iQm,iZc);

    % For each code block
    for cb = 1:maxC
        % Interpolate the code block BLER from the effective SINR using
        % the AWGN table (similar logic to SINR)
        per = interpolatePacketErrorRate(effSINR,iLUT);
    end
end
end

function per = interpolatePacketErrorRate(sinr,lut)
    
end