%% PDSCH Throughput Simulation with CSI Feedback
% First experiment file with CSI-Reporting using RI-PMI-CQI vs Autoencoder
% and calculate throughput
% 
%{

These are the premises of this experiment:
1. Downlink (DL) simulation focus
2. Perfect Uplink (UL)
3. No DL/UL imbalance
4. Uses both Clustered Delay Line (CDL) or Tapped Delay Line (TDL)

%}
% Author: Samnangdona An
% Date: 09/02/2025
% References: 
% 1. https://www.mathworks.com/help/releases/R2025a/5g/ug/nr-pdsch-throughput-using-csi-feedback.html
% 2. https://www.etsi.org/deliver/etsi_ts/138200_138299/138211/16.02.00_60/ts_138211v160200p.pdf
% 3. https://www.sharetechnote.com/html/5G/5G_CSI_RS.html
% 4. https://www.sharetechnote.com/html/5G/5G_CSI_Report.html
% 5. https://www.sharetechnote.com/html/5G/5G_PDSCH_DMRS.html
%
% Data Set provided by The MathWorks, Inc. (www.mathworks.com)
%
% This script and all the functions it used are adapted from reference 1
% files for this specific study, application and setting.

clear
clc
simParams = struct();
simParams.NFrames = 2; % Number of 10ms frames
simParams.SNRIn = -15:5:15;

%% Simulation Toggles
simParams.PerfectChannelEstimator = true;
simParams.EnableCSIRS = false;

%% SCS Settings
simParams.Carrier = nrCarrierConfig;
simParams.Carrier.NSizeGrid = 52;
simParams.Carrier.SubcarrierSpacing = 15;

%% PDSCH
simParams.PDSCH = nrPDSCHConfig;
simParams.PDSCH.PRBSet = 0:simParams.Carrier.NSizeGrid-1;
simParams.PDSCH.NumLayers = 1;
simParams.PDSCH.Modulation = 'QPSK';
simParams.PDSCH.NID = simParams.Carrier.NCellID;

% PDSCH Symbol allocation
% Start with the 3rd index (0-based 2). The first two indices are reserved
% for PDCCH
symAlloc = [2, simParams.Carrier.SymbolsPerSlot-2]; 
simParams.PDSCH.SymbolAllocation = symAlloc;

% DM-RS
%simParams.PDSCH.DMRS.DMRSPortSet = 0:simParams.PDSCH.NumLayers-1;
simParams.PDSCH.DMRS.DMRSTypeAPosition = 2;
simParams.PDSCH.DMRS.DMRSLength = 2;
simParams.PDSCH.DMRS.DMRSAdditionalPosition = 1;
simParams.PDSCH.DMRS.DMRSConfigurationType = 2;
simParams.PDSCH.DMRS.NumCDMGroupsWithoutData = 3;
simParams.PDSCH.DMRS.DMRSEnhancedR18 = false;

% Redundancy Version
simParams.PDSCHExtension.RedundancySequence = 0;

% PDSCH MCS
simParams.PDSCHExtension.MCSTables = nrPDSCHMCSTables;

% LDPC
% Table 5.1.3.1-(1/2/3) TS 38.214
% Picking TargetCodeRate = 490/1024 Intial choice
% This rate is from either Table-1 MCS#13, Table-2 MCS#7, Table-3 MCS#18
MCSIndex = 13;
mcsTable = simParams.PDSCHExtension.MCSTables.QAM64Table;
simParams.LDPC.TargetCodeRate = mcsTable.TargetCodeRate(mcsTable.MCSIndex == MCSIndex); 
simParams.LDPC.DecodingAlgorithm = 'Normalized min-sum';
simParams.LDPC.MaximumLDPCIterationCount = 6;

% Miscellaneous
simParams.PDSCHExtension.PRGBundleSize = 4; % Precoding Block Group
simParams.PDSCHExtension.xOverhead = [];

%% HARQ Process
simParams.EnableHARQ = false;

%% CSI-RS Configuration
simParams.CSIRS = nrCSIRSConfig;
simParams.CSIRS.CSIRSType = 'nzp'; % 'nzp','zp'
simParams.CSIRS.RowNumber = 6; % 1...18
simParams.CSIRS.NumRB = simParams.Carrier.NSizeGrid - simParams.CSIRS.RBOffset;
simParams.CSIRS.CSIRSPeriod = [10 0];
simParams.CSIRS.SymbolLocations = 4;
simParams.CSIRS.SubcarrierLocations = [0,3,6,9];
simParams.CSIRS.Density = 'one';

disp(['Number of CSI-RS ports: ' num2str(simParams.CSIRS.NumCSIRSPorts) '.'])

cdmLengths = getCSIRSCDMLengths(simParams.CSIRS);

%% Channel Model Parameter Settings
simParams.DelayProfile = 'TDL-A';   % 'CDL-' or 'TDL-'
simParams.DelaySpread = 300e-9;     % s
simParams.MaximumDopplerShift = 5;  % Hz

%% Channel Model Compatible Antenna Settings
% Antenna Configuration
% Table of antenna panel array configurations
% M  = # of rows in each antenna panel
% N  = # of columns in each antenna panel
% P  = # of polarizations (1 or 2)
% Mg = # of rows in the array of panels
% Ng = # of columns in the array of panels
% Row format= [N1  N2   P   Mg  Ng]
antArraySizes = ...
   [1   1   1   1   1;   % 1 ants
    1   1   2   1   1;   % 2 ants
    2   1   2   1   1;   % 4 ants
    2   2   2   1   1;   % 8 ants
    2   4   2   1   1;   % 16 ants
    4   4   2   1   1;   % 32 ants
    4   4   2   1   2;   % 64 ants
    4   8   2   1   2;   % 128 ants
    4   8   2   2   2;   % 256 ants
    8   8   2   2   2;   % 512 ants
    8  16   2   2   2];  % 1024 ants

simParams.TransmitAntennaArray.NumPanels        = 1; % Number of transmit panels in horizontal dimension (Ng)
simParams.TransmitAntennaArray.PanelDimensions  = [2, 2]; % Number of columns and rows in the transmit panel (N1, N2)
simParams.TransmitAntennaArray.NumPolarizations = 2; % Number of transmit polarizations
simParams.ReceiveAntennaArray.NumPanels         = 1; % Number of receive panels in horizontal dimension (Ng)
simParams.ReceiveAntennaArray.PanelDimensions   = [2, 1]; % Number of columns and rows in the receive panel (N1, N2)
simParams.ReceiveAntennaArray.NumPolarizations  = 2; % Number of receive polarizations

simParams.NumTxAntennas = prod([simParams.TransmitAntennaArray.NumPanels,...
    simParams.TransmitAntennaArray.PanelDimensions,...
    simParams.TransmitAntennaArray.NumPolarizations]);
simParams.NumRxAntennas = prod([simParams.ReceiveAntennaArray.NumPanels,...
    simParams.ReceiveAntennaArray.PanelDimensions,...
    simParams.ReceiveAntennaArray.NumPolarizations]);

simParams.Channel = Channel.CreateChannel(simParams);
simParams.ChannelInformation = info(simParams.Channel);
fprintf("Simulation Channel Model: %s \n",simParams.DelayProfile);

%% Create Channel clone to calculate Shannon Limit
channelForPathGains = clone(simParams.Channel);
channelForPathGains.ChannelFiltering = false;
channelForPathGains.ChannelResponseOutput = 'path-gains';

%% CSI Report
simParams.CSIReportConfig = struct();
simParams.CSIReportConfig.Mode = 'RI-PMI-CQI';
simParams.CSIReportConfig.Period = [5,0];

switch upper(simParams.CSIReportConfig.Mode)
    case 'RI-PMI-CQI'
        simParams.CSIReportConfig.CQITable          = 'Table1'; % 'Table1','Table2','Table3'
        simParams.CSIReportConfig.CQIMode           = 'Wideband'; % 'Wideband','Subband'
        simParams.CSIReportConfig.PMIMode           = 'Subband'; % 'Wideband','Subband'
        simParams.CSIReportConfig.CodebookType      = 'TypeISinglePanel'; % 'TypeISinglePanel','TypeIIMultiPanel','TypeII','ETypeII'
        simParams.CSIReportConfig.SubbandSize       = 4; % Subband size in RB (4,8,16,32)
        simParams.CSIReportConfig.CodebookMode      = 1; % 1,2
        simParams.CSIReportConfig.RIRestriction     = []; % Empty for no rank restriction
        simParams.CSIReportConfig.NumberOfBeams     = 2; % 2,3,4. Only for Type II codebooks
        simParams.CSIReportConfig.PhaseAlphabetSize = 8; % 4,8. Only for Type II codebooks
        simParams.CSIReportConfig.SubbandAmplitude  = true; % true/false. Only for Type II codebooks
        simParams.CSIReportConfig.ParameterCombination = 1; % 1...8. Only for Enhanced Type II codebooks
        simParams.CSIReportConfig.NumberOfPMISubbandsPerCQISubband = 1; % 1,2. Only for Enhanced Type II codebooks
        simParams.CSIReportConfig.NStartBWP         = []; % Empty to signal the entire carrier
        simParams.CSIReportConfig.NSizeBWP          = []; % Empty to signal the entire carrier
        simParams.CSIReportConfig.RIRestriction = [];

        switch simParams.CSIReportConfig.CodebookType
            case 'TypeIMultiPanel'
                simParams.CSIReportConfig.PanelDimensions = [simParams.ReceiveAntennaArray.NumPanels,...
                    simParams.TransmitAntennaArray.PanelDimensions]; % Panel Dimension
            otherwise
                simParams.CSIReportConfig.PanelDimensions = simParams.TransmitAntennaArray.PanelDimensions; % Panel Dimension
        end
        
    case 'AUTOENCODER'
        aenModelName = 'csiTrainedNetwork.mat';
        if contains(aenModelName,'csiTrainedNetwork')
            simParams.CSIReportConfig.Autoencoder.ModelName = aenModelName;
            simParams.CSIReportConfig.Autoencoder.PretrainedNetwork = true;
            % Load Autoencoder
            [aen,aenOptions,encEndLayer] = Autoencoder.LoadPretrainedNetwork(simParams.CSIReportConfig);
        else
            simParams.CSIReportConfig.Autoencoder.ModelName = aenModelName;
            simParams.CSIReportConfig.Autoencoder.PretrainedNetwork = false;
            % Custom models to be implemented here
            % ...
        end
        % Extract the encoder body
        [csiEncoder,csiDecoder] = Autoencoder.SplitEncoderDecoder(aen,encEndLayer);
        simParams.CSIReportConfig.Autoencoder.Encoder = csiEncoder;
        simParams.CSIReportConfig.Autoencoder.Decoder = csiDecoder;
        simParams.CSIReportConfig.Autoencoder.Options = aenOptions;
    otherwise
        simParams.CSIReportConfig.Mode = 'PERFECT CSI';
end

%% CSI Processsing Delay in Slots
simParams.UEProcessingDelay = 7;
simParams.BSProcessingDelay = 1;

%% Processing Loop
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(simParams.SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(simParams.SNRIn),1);
% Cell array to store CSI reports per SNR point
CSIReportPerSNR = {};
ThisSlotStatus = struct();

for snrIdx = 1:length(simParams.SNRIn)

    % SNR Index
    fprintf('\nSimulating (%dx%d) for SNR = %.2fdB (%s)\n',...
        simParams.NumTxAntennas,simParams.NumRxAntennas,simParams.SNRIn(snrIdx),...
        simParams.DelayProfile);

    % Random seed for reproducibility
    rng(0,'twister');

    % Copy simParams structure to prevent extra variables affecting
    % simParams
    simParamsTmp = simParams;

    % CSI Feedback options
    csiFeedbackOptions = CSIReporting.ConfigureCSIFeedbackOptions(simParamsTmp,snrIdx);

    % Configure Channel
    channel = simParamsTmp.Channel;
    reset(channel);
    reset(channelForPathGains);
    maxChannelDelay = simParamsTmp.ChannelInformation.MaximumChannelDelay;
    nTx = simParamsTmp.ChannelInformation.NumTransmitAntennas;
    nRx = simParamsTmp.ChannelInformation.NumReceiveAntennas;

    % Configure Tx
    [carrier,encDLSCH,pdsch,pdschextra,csirs,wtx] = TxRx.ConfigureTx(simParamsTmp);

    % Configure Rx
    [decDLSCH,timingOffset,N0,noiseEst,csiReports,csiAvailableSlots] = TxRx.ConfigureRx(simParamsTmp,channel,snrIdx,csiFeedbackOptions);

    % Total number of simulation slots
    NSlots = simParamsTmp.NFrames*carrier.SlotsPerFrame;

    % Create an array to store shannon limit
    shannon = zeros(NSlots,1);

    % Loop over waveform length
    for nslot = 0:NSlots-1
        % Update new slot for carrier
        carrier.NSlot = nslot;

        % Determine if there is a new CSI report update
        [isNewCSIReport,reportIdx] = ismember(nslot,csiAvailableSlots);

        % If available, use the new CSI report to configure number of 
        % layers and MCS of the PDSCH
        if isNewCSIReport
            thisCSIReport = csiReports(reportIdx);
            [pdsch.Modulation,pdschextra.TargetCodeRate,wtx] = ...
                CSIReporting.DecodeCSI(carrier,pdsch,pdschextra,thisCSIReport,csiFeedbackOptions);
            pdsch.NumLayers = size(wtx,1);
            encDLSCH.TargetCodeRate = pdschextra.TargetCodeRate;
        end

        % Create resource grid for a slot
        downlinkGrid = nrResourceGrid(carrier,csirs.NumCSIRSPorts);

        % CSI-RS mapping to the slot resource grid
        [csirsIndices,csirsInfo] = nrCSIRSIndices(carrier,csirs);
        csirsSym = nrCSIRS(carrier,csirs);
        downlinkGrid(csirsIndices) = csirsSym;
        isCSIRSOn = ~isempty(csirsIndices);

        % PDSCH reserved REs for CSI-RS
        pdsch.ReservedRE = csirsIndices - 1;

        % Calculate transport block sizes for PDSCH transmission per slot
        [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,...
            length(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,...
            pdschextra.TargetCodeRate,pdschextra.xOverhead);

        % Generate transport block
        for cwIdx = 1:pdsch.NumCodewords
            % Generate new data for current codeword
            x = randi([0,1],trBlkSizes(cwIdx),1);
            encDLSCH.setTransportBlock(x,cwIdx-1);
            decDLSCH.resetSoftBuffer(cwIdx-1);
        end

        % Encode DL-SCH
        %RV = zeros(1,pdsch.NumCodewords);
        xcoded = encDLSCH(pdsch.Modulation,pdsch.NumLayers,...
            pdschIndicesInfo.G,pdschextra.RedundancySequence);

        % Modulate PDSCH and precode
        pdschSyms = nrPDSCH(carrier,pdsch,xcoded);
        % Precode and map to precoding indices
        [pdschPrcSyms,pdschPrcIndices] = nrPDSCHPrecode(carrier,pdschSyms,pdschIndices,wtx);
        downlinkGrid(pdschPrcIndices) = pdschPrcSyms;

        % Add DM-RS
        dmrsSyms = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
        [dmrsPrcSyms,dmrsPrcIndices] = nrPDSCHPrecode(carrier,dmrsSyms,dmrsIndices,wtx);
        downlinkGrid(dmrsPrcIndices) = dmrsPrcSyms;

        % Error if there is any overlapping REs between DM-RS and CSI-RS
        if any(ismember(dmrsIndices,csirsIndices))
            error('CSI-RS and PDSCH DM-RS have overlapping RE(s)');
        end

        % Modulate OFDM Waveform
        txWaveform = nrOFDMModulate(carrier,downlinkGrid);

        % Pass through channel
        % Concatenate with zeros to offset channel delay
        txWaveform = [txWaveform; zeros(maxChannelDelay,nTx)]; %#ok<AGROW>
        [rxWaveform,ofdmResponse,tOffset] = channel(txWaveform,carrier);

        % Get path gain of the channel using the channel clone
        [pathGains,~] = channelForPathGains();
        shannon(nslot+1) = calculateShannonCapacity(carrier,pathGains,nTx,...
            nRx,simParams.SNRIn(snrIdx));

        % Add AWGN to the received waveform
        noisepwr = N0*randn(size(rxWaveform),like=1i);
        rxWaveform = rxWaveform + noisepwr;

        % Assume Perfect timing estimate (no synchronization issue)
        % Timing offset can be added later
        timingOffset = tOffset;
        rxWaveform = rxWaveform(1 + timingOffset:end,:);

        % Demodulate OFDM waveform
        rxDownlinkGrid = nrOFDMDemodulate(carrier,rxWaveform);
        [nRE,nSym,~] = size(rxDownlinkGrid);
        % Pad zero to compensate for channel timing offset shortening of
        % vector
        if nSym < carrier.SymbolsPerSlot
            zeropad = zeros(nRE,carrier.SymbolsPerSlot-nSym,nRx);
            rxDownlinkGrid = cat(2,rxDownlinkGrid,zeropad);
        end

        if simParamsTmp.PerfectChannelEstimator
            % If perfect channel estimator toggle is ON, then the channel
            % response estimate is the ofdmResponse from the channel
            Hest = ofdmResponse;

            % Extract all PDSCH RE from the received downlink grid
            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(...
                pdschIndices,rxDownlinkGrid,Hest);

            % Also precode the channel estimate
            wtxTmp = permute(wtx,[2, 1, 3]); % 2nd-dim is the number of Rx
            pdschHest = nrPDSCHPrecode(carrier,pdschHest,...
                pdschHestIndices,wtxTmp);
        else
            % If perfect channel estimator toggle is OFF, then the channel
            % response estimate is computed via LS algorithm 
            % (nrChannelEstimate) using the DM-RS symbols
            [Hest,noiseEst] = nrChannelEstimate(carrier,rxDownlinkGrid,...
                dmrsIndices,dmrsSyms);

            % Take the average of the noise power across PRGs and layers
            noiseEst = mean(noiseEst,'all');

            % Extract all PDSCH RE from the received downlink grid
            [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxDownlinkGrid,Hest);
        end

        % Equalize the effect of the channel
        [pdschEq,eqCSIScaling] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        % Decode PDSCH
        [dlschLLR,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);

        % Scale LLRs by per-symbol CSI reliability to temper confidence of 
        % noisy subcarriers.
        % Ensures the soft bits fed to the decoder reflect actual SINR post 
        % equalization.
        eqCSIScaling = nrLayerDemap(eqCSIScaling);
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLR{cwIdx})/length(rxSymbols{cwIdx});
            eqCSITmp = repmat(eqCSIScaling{cwIdx},1,Qm).';
            %eqCSIScaling{cwIdx} = repmat(eqCSIScaling{cwIdx},1,Qm);
            dlschLLR{cwIdx} = dlschLLR{cwIdx}.*eqCSITmp(:);
        end

        % Decode DL-SCH
        decDLSCH.TransportBlockLength = trBlkSizes;
        decDLSCH.TargetCodeRate = pdschextra.TargetCodeRate;
        [decBits,blkErr] = decDLSCH(dlschLLR,pdsch.Modulation,...
            pdsch.NumLayers,pdschextra.RedundancySequence);

        % Save values to calculate throughput
        simThroughput(snrIdx) = simThroughput(snrIdx) + sum(~blkErr .* trBlkSizes);
        maxThroughput(snrIdx) = maxThroughput(snrIdx) + sum(trBlkSizes);

        % CSI measurement and encoding
        if isCSIRSOn

            % Generate CSI report
            rxCSIReport = CSIReporting.EncodeCSI(carrier,csirs,Hest,noiseEst,csiFeedbackOptions);
            reportPeriod = csiFeedbackOptions.CSIReportConfig.Period(1);
            reportOffset = csiFeedbackOptions.CSIReportConfig.Period(2);
            reportNSlot  = 1 + nslot + simParamsTmp.UEProcessingDelay; % Accounts for UE proc delay
            csiFeedbackSlot = stepNextCSISlot(reportPeriod,reportOffset,reportNSlot);
            csiAvailableSlots(end+1) = 1 + csiFeedbackSlot + simParamsTmp.BSProcessingDelay; %#ok<SAGROW> Accounts for BS proc delay
            csiReports(end+1) = rxCSIReport; %#ok<SAGROW>

        end
        % Aggregate slot transmission status
        ThisSlotStatus.SlotNumber = nslot;
        ThisSlotStatus.HasBlockError = blkErr;
        ThisSlotStatus.PDSCH = pdsch;
        ThisSlotStatus.PDSCHExtra = pdschextra;
        ThisSlotStatus.CodeRate = trBlkSizes./pdschIndicesInfo.G;
        ThisSlotStatus.IsCSIRSOn = isCSIRSOn;
        ThisSlotStatus.CSIReportIndex = reportIdx;
        ThisSlotStatus.CSIReport = csiReports;

        % Display slot status info
        displaySlotTransmissionStatus(ThisSlotStatus,NSlots,carrier);

    end

    % Store CSI Report for each SNR
    CSIReportPerSNR{snrIdx} = csiReports;  %#ok<SAGROW>
    
    % Display throughput for each SNR
    displayThroughput(simParamsTmp,snrIdx,simThroughput,shannon);

end

%% Local Functions
function cdmLengths = getCSIRSCDMLengths(csirs)
cdm = csirs.CDMType;
opt = {  'noCDM',[1,1]
       'fd-CDM2',[2,1]
          'CMD4',[2,2]
          'CDM8',[2,4]};
cdmLengths = opt{strcmpi(cdm,opt(:,1)),2};
end

function csiSlot = stepNextCSISlot(period,offset,nSlot)
% Output the slot number of the next suitable slot for CSI reporting per
% CSI report config periodicity

csiSlot = period*ceil((nSlot-offset)/period)+offset;
end

function [C,SE] = calculateShannonCapacity(carrier,HGrid,nTx,nRx,SNRdB)
%calculateShannonCapacityUpperBound calculates the Shannon limit for the
%   MIMO communication link.
%
%   C = B * log2(det(I_nRx + (SNR/nTx) * H * H_hermitian))
%
% Inputs:
%   carrier : struct/object with NSizeGrid and SubcarrierSpacing (kHz)
%   H       : nRx x nTx channel matrix (single realization)
%   nTx     : number of Tx antennas
%   nRx     : number of Rx antennas
%   SNR_dB  : SNR in dB (per receive antenna)

% Bandwidth = N_PRB * NumSubcarrierPerPRB * SubcarrierSpacing (Hz)
B = carrier.NSizeGrid * 12 * carrier.SubcarrierSpacing * 1e3;

% Get effective channel response
HavgTime = squeeze(mean(HGrid,1)); % [PathDelays x nTx x nRx]
H = squeeze(sum(HavgTime,1)); % sum over delays -> [nTx x nRx]
H = permute(H,[2, 1]); % Make it [nRx x nTx]
HHt = H*H'; % H * H^H (nRx x nTx)

% Eigenvalue-based capacity (more stable than determinant)
lambda = eig(HHt);  % nRx x 1
SNR = 10.^(SNRdB/10);
SE = sum(log2(1 + (SNR/nTx)*lambda)); % spectral efficiency
C = B*SE;

% Alternatively (strictly from formula)
% InRx = eye(nRx);
% M = InRx + (SNR/nTx)*HHt;
% C = B * log2(real(det(M)));
end

function displaySlotTransmissionStatus(slotStatus,TotalNumSlots,carrier)

pdsch = slotStatus.PDSCH;
pdschextra = slotStatus.PDSCHExtra;
thisSlotCSIReport = slotStatus.CSIReport;
thisSlotCSIReportIndex = slotStatus.CSIReportIndex;
numCodewords = pdsch.NumCodewords;
isMultiCodewords = numCodewords > 1;
codewordLayers = floor((pdsch.NumLayers + (0:numCodewords-1))/numCodewords);
tcr = pdschextra.TargetCodeRate; % target code rate
acr = slotStatus.CodeRate; % actual code rate
modulation = pdsch.Modulation;

for cwIdx = 1:numCodewords
    if slotStatus.HasBlockError
        txStatus = 'Transmission Failed';
    else
        txStatus = 'Transmission Succeeded';
    end
    infoStringCW = sprintf('%s (Layers=%d, %s, TargetCR=%.3f, ActualCR=%.3f)',...
        txStatus,codewordLayers(cwIdx),modulation{cwIdx},tcr(cwIdx),acr);

    if isMultiCodewords
        infoToDisplay = sprintf('Codeword #%d %s',cwIdx-1,infoStringCW);
    else
        infoToDisplay = infoStringCW;
    end
end

csirsInfoStr = [];
if slotStatus.IsCSIRSOn
    csirsInfoStr = "CSI-RS Transmission Active. ";
end

csiFeedBackInfoStr = [];
nslot = carrier.NSlot;
if nslot == 0
    % First ever slot
    csiFeedBackInfoStr = "Using initial CSI | ";
elseif thisSlotCSIReportIndex > 0
    csiFeedBackInfoStr = sprintf("Using CSI from NSlot: #%2d",thisSlotCSIReport(thisSlotCSIReportIndex).NSlot);
end

strtmp = join([infoToDisplay,csiFeedBackInfoStr,csirsInfoStr]);
fprintf("(%5.2f%%) NSlot: #%2d: %s \n",100*(nslot+1)/TotalNumSlots,nslot,strtmp);

end

function displayThroughput(simParams,snrIdx,throughput,shannon)

frames = simParams.NFrames;
snrVal = simParams.SNRIn(snrIdx);
mbps = (throughput(snrIdx)/(simParams.NFrames*10e-3))/1e6;

fprintf('\nThroughput at SNR = %.2fdB for %d frames: %.3f Mbps\n',snrVal,frames,mbps);
fprintf('The average Shannon capacity (ergodic Shannon limit) of this link: %.3f Mbps\n',mean(shannon)/1e6);

end


