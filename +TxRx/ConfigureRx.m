function [decDLSCH,timingOffset,N0,noiseEst,csiReports,csiAvailableSlots] = ...
    ConfigureRx(simParams,channel,snrIdx,csiFeedbackOpts)
% Create and configure DL-SCH decoder. Obtain noise related quantities and
% initial CSI feedback from perfect channel knowledge.

decDLSCH = nrDLSCHDecoder;
decDLSCH.LDPCDecodingAlgorithm = simParams.LDPC.DecodingAlgorithm;
decDLSCH.MaximumLDPCIterationCount = simParams.LDPC.MaximumLDPCIterationCount;
decDLSCH.MultipleHARQProcesses = simParams.EnableHARQ;

% Carrier and waveform information
carrier  = simParams.Carrier;
slotsPerSubFrame = carrier.SlotsPerSubframe;
waveformInfo = nrOFDMInfo(carrier);
nfft = waveformInfo.Nfft;
sampleRate = waveformInfo.SampleRate;

% --- Noise power estimation block ---
% Goal: compute noise variance consistent with target SNR, FFT size, and
% channel scaling (especially for multi-antenna configs).
% Steps:
% 1. Convert desired SNR from dB to linear.
% 2. Derive base noise std-dev (N0) per FFT bin.
%    - Scaled by sqrt(nfft) to normalize for OFDM FFT/IFFT energy spread.
%    - Scaled by sqrt(SNRLin) so that resulting signal/noise ratio matches
%      the requested target.
% 3. Adjust for multiple Tx antennas if channel outputs are normalized.
% 4. Convert noise std-dev into noise *power* across the FFT grid.
SNRdB  = simParams.SNRIn(snrIdx);       % Desired SNR [dB]
SNRLin = 10^(SNRdB/10);                 % Convert to linear
N0 = 1/sqrt(nfft*SNRLin);               % Base noise scaling (per RE, per FFT bin)

numInputs = prod(channel.TransmitAntennaArray.Size);  % Total Tx antenna elems
if channel.NormalizeChannelOutputs
    N0 = N0/sqrt(numInputs);            % Compensate for normalized channel outputs
end
noiseEst = N0^2*nfft;                   % Final noise variance (per FFT grid)

% Channel Info
channelInfo = info(channel);

% Get initial channel estimate
channelClone = clone(channel); % clone to get ofdm-response without affect original object
release(channelClone)
channelClone.ChannelFiltering = false;
channelClone.ChannelResponseOutput = 'ofdm-response';
maxChannelDelay = channelInfo.MaximumChannelDelay; % ceil(max(pathDelays*sampleRate)) + filterDelay
channelClone.NumTimeSamples = (sampleRate/1000)/slotsPerSubFrame + maxChannelDelay;
[Hest,timingOffset] = channelClone(carrier); % Hest is the ofdm-response

% Get CSI-RS initial report based on perfect channel estimate
csirs = simParams.CSIRS;
csiFeedbackOpts.PerfectChannelEstimator = true;
csirs.CSIRSPeriod = 'on';
csiReports = CSIReporting.EncodeCSI(carrier,csirs,Hest,noiseEst,csiFeedbackOpts);
csiAvailableSlots = 0;

end