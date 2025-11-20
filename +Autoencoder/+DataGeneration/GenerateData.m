function [HTargetTensor,HTargetTensorpp,options] = GenerateData(numSamples,carrier,channel,options)
%GenerateData produces N channel estimate frames using the 5G NR carrier
%   (CR) and channel model (CH). The output HTarget is a tensor of shape 
%   [Ndelay × NTx × 2(I/Q) × Nrx], where:
%       - Ndelay = max channel delay taps
%       - Ntx = transmit antennas
%       - 2 = I/Q components
%       - Nrx = receive antennas
% opt holds configuration parameters. 
% Each saved file corresponds to one channel estimate frame.

release(channel)
channelInfo = info(channel);
switch class(channel)
    case 'nrCDLChannel'
        nTx = channelInfo.NumInputSignals;
        nRx = channelInfo.NumOutputSignals;
        options.ChannelSampleDensity = 64*4;
    case 'nrTDLChannel'
        nTx = channelInfo.NumInputSignals;
        nRx = channelInfo.NumReceiveAntennas;
        options.ChannelSampleDensity = 64;
    otherwise
        error('Unknown or unsupported channel type')
end

% Waveform information
waveInfo = nrOFDMInfo(carrier);
channel.SampleRate = waveInfo.SampleRate; % set same sample rate

% Samples per slot based on carrier waveform 
symbolsPerSlot = carrier.SymbolsPerSlot;
options.SamplesPerSlot = sum(waveInfo.SymbolLengths(1:symbolsPerSlot));

% Total number of subcarriers in the entire resource grid
numSubcarriers = carrier.NSizeGrid*12; % 12 subcarrier per RB

% Frequency resolution per subcarrier
freqResolution = carrier.SubcarrierSpacing*1000; % 1000 to convert kHz->Hz

% Total occupied bandwidth per resource grid
totalBandwidth = numSubcarriers*freqResolution; % in Hz

% Delay or Time resolution: the smallest distinguishable difference between
% two multipath delays
delayResolution = 1/totalBandwidth;

% Max delay
options.MaxDelay = round((channel.DelaySpread/delayResolution)*options.TruncationFactor/2)*2;

% Various opt params
options.NumSlotsPerFrame = 1;
options.Preprocess = true;
options.ResetChannelPerFrame = true;
options.Normalization = false;

% Calculate the number of required frames to accomodate data generation for
% the given configurations
% Each slot has nRx number of training channel estimates.
numSlots = ceil(numSamples/nRx);
numFrames = ceil(numSlots/options.NumSlotsPerFrame);

if options.SaveData
    fileNamePrefix = opt.DataFilePrefix;
    HTargetFileNameBase = fullfile(pwd,opt.DataDir,fileNamePrefix);
    if ~exist(fullfile(pwd,opt.DataDir),"file")
        mkdir(fullfile(pwd,opt.DataDir))
    end
    if opt.Preprocess
        processedFileNameBase = ...
            fullfile(pwd,opt.DataDir,"processed",fileNamePrefix+"_processed");
        if ~exist(fullfile(pwd,opt.DataDir,"processed"),"file")
            mkdir(fullfile(pwd,opt.DataDir,"processed"))
        end
    else
        processedFileNameBase = [];
    end
else
    HTargetFileNameBase = [];
    processedFileNameBase = [];
end

% Pre-allocate target tensors
HTargetTensor = zeros(numSubcarriers,options.NumSlotsPerFrame*symbolsPerSlot,nRx,nTx,numFrames,like=1i);
HTargetTensorpp = zeros(options.MaxDelay,nTx,2,nRx,options.NumSlotsPerFrame*numFrames);

for frIdx = 1:numFrames

    if options.Verbose && (rem(frIdx,100) == 0 || frIdx == 1 || frIdx == numFrames)
        fprintf('Processing frame #%d of %d \n',frIdx,numFrames);
    end

    channelTmp = clone(channel);

    % Generate channel response tensor
    HestTensor = estimateChannelResponse(channelTmp,carrier,options);

    % Fill in H target tensor
    HTargetTensor(:,:,:,:,frIdx) = HestTensor;

    %save(processedFileNameBase + "_" + fr,'HestTensor');

    % Preprocess channel response tensor
    HestTensorpp = preprocessChannelResponse(HestTensor,options);

    % Fill in preprocessed H target tensor
    HTargetTensorpp(:,:,:,:,frIdx) = HestTensorpp;
end
fprintf('Data Generation Completed \n')

end


%% Local Helper fcn
function HestTensor = estimateChannelResponse(channel,carrier,options)
%EstimateChannelResponse Channel response estimation for data generation.
%   Generates one frame of perfect channel estimate that contains N slots
%   per frame.

% Perform slot and frame compatibility check
% Is the number of slots we want to simulate smaller than one full 5G 
% frame?
slotsPerFrame = options.NumSlotsPerFrame;
symsPerSlot = carrier.SymbolsPerSlot;
channel.SampleDensity = options.ChannelSampleDensity;

if slotsPerFrame < symsPerSlot
    % If smaller, we're only simulating less than one full frame
    channel.NumTimeSamples = options.SamplesPerSlot*slotsPerFrame;

    % In this case, only generate once from the subroutine
    HestTensor = channelEstimationSubroutine(channel,carrier,options.ZeroTimingOffset);
else
    % If greater, we'll have to repeatedly generate channel estimation and
    % stitch them together

    % Calculate the ratio of slotsPerFrame over carrier's slotsPerFrame
    ratio = slotsPerFrame/carrier.SlotsPerFrame;

    % Calculate the total number of symbols
    totalNumSyms = symsPerSlot*ceil(ratio)*carrier.SlotsPerFrame;

    % Pre-allocate Hest tensor
    HestTensor = zeros(options.NumSubcarriers,totalNumSyms,options.nRx,options.nTx,like=1i);

    for subFrameIdx = 0:ceil(ratio)-1
        HestTmp = channelEstimationSubroutine(channel,carrier,options.ZeroTimingOffset);
        symbolsPerFrame = carrier.SlotsPerFrame*carrier.SymbolsPerSlot;
        startIdx = subFrame*symbolsPerFrame + 1;
        stopIdx  = (subFrame + 1)*symbolsPerFrame;
        HestIdx  = startIdx:stopIdx;

        % Load HestTmp into Hest tensor
        HestTensor(:,HestIdx,:,:) = HestTmp;
    end

    if opt.ResetChannelPerFrame
        reset(channel)
    end
end
end

function Hest = channelEstimationSubroutine(channel,carrier,zto)
%channelEstimationSubroutine subroutine to help with channele response
%   estimation for different slots

[pathGains,sampleTimes] = channel();
pathFilters = channel.getPathFilters();

if zto
  % Perfect timing sync
  offset = 0;
else
  offset = nrPerfectTimingEstimate(pathGains,pathFilters);
end

% Have to use perfect channel estimate for targets/ground truth
Hest = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

end

function Hestpp = preprocessChannelResponse(Hest,options)
%preprocessChannelResponse Prepares the channel estimate 'Hest' by
%   decimating via 2D-FFT and 2D-IFFT. The preprocessed channel estimate
%   has maxDelay samples and is a real-valued tensor.
%   Hestpp is maxDelay x nTx x 2(I/Q) x nRx x Nslots, where 
%   Nslots = numSyms/14
%
%   This procedure is similar to Principle Components Analysis (PCA)
%

maxDelay = options.MaxDelay;
[numSubcarrier,numSyms,nRx,nTx] = size(Hest);
% Compute center of delay spread band. Keep only the taps that matter and
% discard channel where magnitude is close to 0. This is delay-domain
% sparsity exploitation.
midpt = floor(numSubcarrier/2);
lb = midpt - (numSubcarrier - maxDelay)/2 + 1;
ub = midpt + (numSubcarrier - maxDelay)/2;

% Average over symbols
symsPerSlot = 14;
Nslots = numSyms/symsPerSlot;
HTmp = reshape(Hest,numSubcarrier,Nslots,symsPerSlot,nRx,nTx);
% Remove fast fading noise and keeps long-term multipath structure
H = mean(HTmp,3); % 3rd dim is symsPerSlot dim
H = permute(H, [1, 5, 4, 2, 3]); % rearrange dims

% Pre-allocate
Hestpp = zeros(maxDelay,nTx,2,nRx,Nslots,like=1i);

% Per Rx antenna, decimate over subcarriers using 2D-FFT
for slotIdx = 1:Nslots
    % 2D-FFT: go from freq/spatial(antenna) to delay/angle(spatial dir) 
    % domain
    Hfft = fft2(H(:,:,:,slotIdx));

    % Truncate
    HTruncTmp = Hfft([1:lb-1, ub+1:end],:,:);

    switch options.DataDomain
        case 'Frequency-Spatial'
            % Transform back if freq/spatial is desired
            HTrunc = ifft2(HTruncTmp);
        otherwise
            % Keep as is iff delay/spatial is desired
            HTrunc = HTruncTmp;
    end

    % Extract Real (I) and Imag (Q) components
    HI = real(HTrunc); % I
    HQ = imag(HTrunc); % Q

    if options.Normalization
        avgVal = options.MeanValue;
        stdVal = options.StdValue;
        targetStd = options.TargetStdValue;
        % Normalize using standard Z-score and rescaling std
        Hestpp(:,:,1,:,slotIdx) = ((HI - avgVal)/stdVal)*targetStd + 0.5;
        Hestpp(:,:,2,:,slotIdx) = ((HQ - avgVal)/stdVal)*targetStd + 0.5;
    else
        % Keep as is
        Hestpp(:,:,1,:,slotIdx) = HI;
        Hestpp(:,:,2,:,slotIdx) = HQ;
    end
end

end
















