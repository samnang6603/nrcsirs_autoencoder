function channel = createChannel(simParams)
%CREATECHANNEL Create propagation channel based on simParams. 
% Channel can be CDL or TDL

% Get waveform information, specifically the sample rate
waveInfo = nrOFDMInfo(simParams.Carrier);
sampleRate = waveInfo.SampleRate;

% Get common propagation channel parameters
delayProfile = simParams.DelayProfile;
delaySpread = simParams.DelaySpread;
maximumDopplerShift = simParams.MaximumDopplerShift;

% Number of antenna elements and polarizations
if contains(delayProfile,'CDL')
    % Construct CDL channel 
    channel = nrCDLChannel(DelayProfile = delayProfile,...
        DelaySpread = delaySpread,...
        MaximumDopplerShift = maximumDopplerShift,...
        ChannelResponseOutput = 'ofdm-response',...
        SampleRate = sampleRate);

    % CDL requires configurationss of antenna structure
    % Polarization is always 2 if number of antenna element > 2
    nTxPol = (simParams.NTxAnts > 1) + 1;
    nRxPol = (simParams.NRxAnts > 1) + 1;

    % Tx antenna array configuration
    txArray = simParams.TransmitAntennaArray;
    M = txArray.PanelDimensions(2);
    N = txArray.PanelDimensions(1);
    Ng = txArray.NumPanels;

    channel.TransmitAntennaArray.Size = [M N nTxPol 1 Ng];
    channel.TransmitAntennaArray.ElementSpacing = [0.5 0.5 1 1]; % Element spacing in wavelengths
    channel.TransmitAntennaArray.PolarizationAngles = [-45 45];  % Polarization angles in degrees

    % Rx antenna array configuration
    rxArray = simParams.ReceiveAntennaArray;
    M = rxArray.PanelDimensions(2);
    N = rxArray.PanelDimensions(1);
    Ng = rxArray.NumPanels;

    channel.ReceiveAntennaArray.Size = [M N nRxPol 1 Ng];
    channel.ReceiveAntennaArray.ElementSpacing = [0.5 0.5 1 1];  % Element spacing in wavelengths
    channel.ReceiveAntennaArray.PolarizationAngles = [0 90];     % Polarization angles in degrees
elseif contains(delayProfile,'TDL')
    % Construct TDL Channel
    channel = nrTDLChannel(DelayProfile = delayProfile,...
        DelaySpread = delaySpread,...
        MaximumDopplerShift = maximumDopplerShift,...
        ChannelResponseOutput = 'ofdm-response',...
        SampleRate = sampleRate);

    % TDL only requires total number of Tx and Rx
    channel.NumTransmitAntennas = simParams.NTxAnts;
    channel.NumReceiveAntennas  = simParams.NRxAnts;
end