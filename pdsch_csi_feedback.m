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
% This script and all the functions it used are adapted from reference 1
% files for this specific study, application and setting.

clear
clc
simParams = struct();
simParams.SNRIn = -10:5:10;

%% Simulation Toggles
simParams.PerfectChannelEstimator = false;
simParams.EnableCSIRS = false;

%% SCS Settings
simParams.Carrier = nrCarrierConfig;
simParams.Carrier.NSizeGrid = 52;
simParams.Carrier.SubcarrierSpacing = 15;

%% PDSCH
simParams.PDSCH = nrPDSCHConfig;
simParams.PDSCH.PRBSet = 0:simParams.Carrier.NSizeGrid-1;
simParams.PDSCH.NumLayers = 2;
simParams.PDSCH.Modulation = 'QPSK';
simParams.PDSCH.NID = simParams.Carrier.NCellID;

% DM-RS
simParams.PDSCH.DMRS.DMRSPortSet = 0:simParams.PDSCH.NumLayers-1;
simParams.PDSCH.DMRS.DMRSTypeAPosition = 2;
simParams.PDSCH.DMRS.DMRSLength = 2;
simParams.PDSCH.DMRS.DMRSAdditionalPosition = 1;
simParams.PDSCH.DMRS.DMRSConfigurationType = 1;
simParams.PDSCH.DMRS.NumCDMGroupsWithoutData = 2;

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

%% Antenna Configuration
% Table of antenna panel array configurations
% M  = # of rows in each antenna panel
% N  = # of columns in each antenna panel
% P  = # of polarizations (1 or 2)
% Mg = # of rows in the array of panels
% Ng = # of columns in the array of panels
% Row format= [M  N   P   Mg  Ng]
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
simParams.TransmitAntennaArray.PanelDimensions  = [1, 1]; % Number of columns and rows in the transmit panel (N1, N2)
simParams.TransmitAntennaArray.NumPolarizations = 2; % Number of transmit polarizations
simParams.ReceiveAntennaArray.NumPanels         = 1; % Number of receive panels in horizontal dimension (Ng)
simParams.ReceiveAntennaArray.PanelDimensions   = [1, 1]; % Number of columns and rows in the receive panel (N1, N2)
simParams.ReceiveAntennaArray.NumPolarizations  = 2; % Number of receive polarizations

simParams.NTxAnts = prod([simParams.TransmitAntennaArray.NumPanels,...
    simParams.TransmitAntennaArray.PanelDimensions,...
    simParams.TransmitAntennaArray.NumPolarizations]);
simParams.NRxAnts = prod([simParams.ReceiveAntennaArray.NumPanels,...
    simParams.ReceiveAntennaArray.PanelDimensions,...
    simParams.ReceiveAntennaArray.NumPolarizations]);

%% CSI-RS Configuration
simParams.CSIRS = nrCSIRSConfig;
simParams.CSIRS.CSIRSType = 'nzp'; % 'nzp','zp'
simParams.CSIRS.RowNumber = 3; % 1...18
simParams.CSIRS.NumRB = simParams.Carrier.NSizeGrid - simParams.CSIRS.RBOffset;
simParams.CSIRS.CSIRSPeriod = [10 0];
simParams.CSIRS.SymbolLocations = 4;
simParams.CSIRS.SubcarrierLocations = 0; %[0,3,6,9];
simParams.CSIRS.Density = 'one';

disp(['Number of CSI-RS ports: ' num2str(simParams.CSIRS.NumCSIRSPorts) '.'])

cdmLengths = getCSIRSCDMLengths(simParams.CSIRS);

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
        simParams.Autoencoder = 'csiTrainedNetwork.mat';
    otherwise
        simParams.CSIReportConfig.Mode = 'PERFECT CSI';
end

%% CSI Processsing Delay in Slots
simParams.UEProcessingDelay = 7;
simParams.BSProcessingDelay = 1;

%% Channel
simParams.DelayProfile = 'CDL-C';   % 'CDL-' or 'TDL-'
simParams.DelaySpread = 300e-9;     % s
simParams.MaximumDopplerShift = 5;  % Hz
simParams.Channel = Channel.CreateChannel(simParams);
simParams.ChannelInformation = info(simParams.Channel);

%% Processing Loop
% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(simParams.SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(simParams.SNRIn),1);
% Cell array to store CSI reports per SNR point
CSIReport = {};

for snrIdx = 1:length(simParams.SNRIn)

    % Random seed for reproducibility
    rng(0,'twister');

    % Copy simParams structure to prevent extra variables affecting
    % simParams
    simParamsPrivate = simParams;

    % CSI Feedback options
    csiFeedbackOptions = CSIReporting.ConfigureCSIFeedbackOptions(simParamsPrivate,snrIdx);

    % Configure Channel
    channel = simParamsPrivate.Channel;
    reset(channel);
    maxChannelDelay = simParamsPrivate.ChannelInformation.MaximumChannelDelay;

    % Configure Tx
    [carrier,encDLSCH,pdsch,pdschextra,csirs,wtx] = TxRx.ConfigureTx(simParamsPrivate);

    % Configure Rx
    [decodeDLSCH,timingOffset,N0,noiseEst,csiReports,csiAvailableSlots] = TxRx.ConfigureRx(simParamsPrivate,channel,snrIdx,csiFeedbackOptions);

    

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





