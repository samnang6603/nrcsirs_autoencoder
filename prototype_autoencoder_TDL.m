% Prototype of Autoencoder for CSI-RS Report Feedback
% Reference:
% 1. https://www.mathworks.com/help/releases/R2025a/5g/ug/csi-feedback-with-autoencoders.html
% 2. Zimaglia, Elisa, Daniel G. Riviello, Roberto Garello, and Roberto Fantini.
%   "A Novel Deep Learning Approach to CSI Feedback Reporting for NR 5G Cellular Systems."
%   In 2020 IEEE Microwave Theory and Techniques in Wireless Communications (MTTW),
%   47–52. Riga, Latvia: IEEE, 2020.
%   https://doi.org/10.1109/MTTW51045.2020.9245055.

clear, clc

%% Architecture

trainingToggle = true;

encoderLayers = [
    imageInputLayer([28, 2, 2],Normalization='none',Name='Enc_Input')
    convolution2dLayer([3, 3],2,Padding='same',Name='Enc_Conv')
    batchNormalizationLayer(Epsilon=0.001,Name='Enc_BN')
    leakyReluLayer(0.3,Name='Enc_LeakyRelu')
    flattenLayer(Name='Enc_Flatten')
    fullyConnectedLayer(64,Name='Enc_FCN')
    sigmoidLayer(Name='Enc_Sig')
    ];

tmp = [fullyConnectedLayer(28*2*2,Name='Dec_FCN');
    functionLayer(@(x)dlarray(reshape(x,28,2,2,[]),'SSCB'),Formattable=true,...
    Acceleratable=true,Name='Dec_Reshape_S1')];

aen = dlnetwork([encoderLayers; tmp]);

decoderStage1Layers = [
    convolution2dLayer([3, 3],8,Padding='same',Name='Dec_Conv_S1_1')
    batchNormalizationLayer(Epsilon=0.001,Name='Dec_BN_S1_1')
    leakyReluLayer(0.3,Name='Dec_LeakyRelu_S1_1')

    convolution2dLayer([3, 3],8*2,Padding='same',Name='Dec_Conv_S1_2')
    batchNormalizationLayer(Epsilon=0.001,Name='Dec_BN_S1_2')
    leakyReluLayer(0.3,Name='Dec_LeakyRelu_S1_2')

    convolution2dLayer([3 3],2,Padding='same',Name='Dec_Conv_S1_3')
    batchNormalizationLayer(Epsilon=0.001,Name='BN_S1_3')

    ];

aen = addLayers(aen,decoderStage1Layers);
aen = connectLayers(aen,'Dec_Reshape_S1','Dec_Conv_S1_1');
addStage1Layer = additionLayer(2,Name='Add_S1_oBN3_To_S1_oReshape');
aen = addLayers(aen,addStage1Layer);
aen = connectLayers(aen,'BN_S1_3','Add_S1_oBN3_To_S1_oReshape/in1');
aen = connectLayers(aen,'Dec_Reshape_S1','Add_S1_oBN3_To_S1_oReshape/in2');
finalLeakyReluStage1Layer = leakyReluLayer(0.3,Name='Dec_LeakyRelu_S1_3');
aen = addLayers(aen,finalLeakyReluStage1Layer);
aen = connectLayers(aen,'Add_S1_oBN3_To_S1_oReshape','Dec_LeakyRelu_S1_3');

decoderStage2Layers = [

convolution2dLayer([3, 3],8,Padding='same',Name='Dec_Conv_S2_1')
batchNormalizationLayer(Epsilon=0.001,Name='Dec_BN_S2_1')
leakyReluLayer(0.3,Name='Dec_LeakyRelu_S2_1')

convolution2dLayer([3, 3],8*2,Padding='same',Name='Dec_Conv_S2_2')
batchNormalizationLayer(Epsilon=0.001,Name='Dec_BN_S2_2')
leakyReluLayer(0.3,Name='Dec_LeakyRelu_S2_2')

convolution2dLayer([3 3],2,Padding='same',Name='Dec_Conv_S2_3')
batchNormalizationLayer(Epsilon=0.001,Name='BN_S2_3')

];

aen = addLayers(aen,decoderStage2Layers);
aen = connectLayers(aen,'Dec_LeakyRelu_S1_3','Dec_Conv_S2_1');
addStage2Layer = additionLayer(2,Name='Add_S1_oLeakyRelu3_To_S2_oBN3');
aen = addLayers(aen,addStage2Layer);
aen = connectLayers(aen,'BN_S2_3','Add_S1_oLeakyRelu3_To_S2_oBN3/in1');
aen = connectLayers(aen,'Dec_LeakyRelu_S1_3','Add_S1_oLeakyRelu3_To_S2_oBN3/in2');
finalLeakyReluStage2Layer = leakyReluLayer(0.3,Name='Dec_LeakyRelu_S2_3');
aen = addLayers(aen,finalLeakyReluStage2Layer);
aen = connectLayers(aen,'Add_S1_oLeakyRelu3_To_S2_oBN3','Dec_LeakyRelu_S2_3');

finalConvSig = [
    convolution2dLayer([3 3],2,Padding='same',Name='Dec_Conv_Final')
    sigmoidLayer(Name='Dec_Sig_Final')
    ];

aen = addLayers(aen,finalConvSig);
aen = connectLayers(aen,'Dec_LeakyRelu_S2_3','Dec_Conv_Final');

%figure
%plot(aen)
%title('Autoencoder Prototype Structure (MATLAB Example)')

%% Carrier
simParams.Carrier = nrCarrierConfig;

%% Create Channel
channel = nrTDLChannel;
channel.DelayProfile = 'TDL-A';
channel.DelaySpread = 300e-9;       % s
channel.MaximumDopplerShift = 5;   % Hz
channel.RandomStream = "Global stream";
channel.NumTransmitAntennas = 2;
channel.NumTransmitAntennas = 2;
channel.ChannelFiltering = false;           % No filtering

fprintf("Generate %s channel response \n",channel.DelayProfile);

%% Autoencoder Options
aenOptions.DataDomain = 'Frequency-Spatial';
aenOptions.TruncationFactor = 10;
aenOptions.ZeroTimingOffset = true;
aenOptions.SaveData = false;
aenOptions.Verbose = true;

%% Generate channel response data
%rng(0)
numSamples = 1500;
[~,Hpp,aenOptions] = Autoencoder.DataGeneration.GenerateData( ...
    numSamples,simParams.Carrier,channel,aenOptions);

% Each frame has a data for nRx
[maxDelay,nTx,NIQ,nRx,Nframes] = size(Hpp);

% Combine frames and nRx
Hpp = reshape(Hpp,maxDelay,nTx,NIQ,nRx*Nframes);

% Calculate the average and standard deviation to normalize the data
avgVal = mean(Hpp,'all');
stdVal = std(Hpp,[],'all');

% Separate data into training, validation and testing
numDataSet = size(Hpp,4);
numTrainDS = (10/15)*numDataSet; %0.6*numDataSet;  % 60%
numValidDS = (3/15)*numDataSet; %0.15*numDataSet; % 15%
numTestDS  = (2/15)*numDataSet; %0.25*numDataSet; % 25%

% Target Std of value of 0.0212, which will restrict data to range
% [-0.5 , 0.5]
trainInd = 1:numTrainDS;
validInd = numTrainDS + (1:numValidDS);
testInd  = numTrainDS + numValidDS + (1:numTestDS);
HTrain = (Hpp(:,:,:,trainInd) - avgVal)/stdVal*0.0212 + 0.5;
HValid = (Hpp(:,:,:,validInd) - avgVal)/stdVal*0.0212 + 0.5;
HTest  = (Hpp(:,:,:,testInd)  - avgVal)/stdVal*0.0212 + 0.5;
if ispc
    if ~isfolder("nrcsirs_autoencoder\Data")
        mkdir("Data\TDL\");
    end
    save("Data\TDL\" + "ChannelEstimationTDL.mat",...
        'Hpp','HTrain','HTest','HValid');
end

% Update aen options
aenOptions.AverageValue = avgVal;
aenOptions.StandardDeviationValue = stdVal;
aenOptions.TargetStandardDeviationValue = 0.0212;

%% Training
if trainingToggle
    % Training options
    miniBatch = 1000;
    valData = {HValid,HValid};
    [trainingOpts,lossFcn] = Autoencoder.Training.ConfigureOptions('adam',valData,miniBatch);

    % The training data doubles as the ground truth, since the autoencoder
    % learns to compress the input and then reconstruct that same input at
    % the output.
    HTarget = HTrain;

    % Start training
    disp('Begin Training ... ')
    [trainedNet,trainInfo] = trainnet(HTrain,HTarget,aen,lossFcn,trainingOpts);
    savedTrainingOpts = trainingOpts;
    if ispc
        if ~isfolder("nrcsirs_autoencoder\Saved_Models")
            mkdir("Saved_Models\TDL\");
        end
        ss = "AENTrainedNetwork_" + ...
            string(datetime('now',Format="MM_dd_HH_mm"));
        save("Saved_Models\TDL\" + ss,...
            'trainedNet','trainInfo','aenOptions','savedTrainingOpts');
        fprintf("Model %s has been saved to Saved_Models folder \n",ss);
    end

else
    [fname, fpath, ~] = uigetfile("*.mat",MultiSelect="on");
    load([fpath,fame]);
end


%% Test dataset
HTestHat = predict(trainedNet,HTest);

% Get performance stats
correlation = zeros(numTestDS,1);
nmse = zeros(numTestDS,1);
for idx = 1:numTestDS
    x = HTest(:,:,1,idx) + 1i*HTest(:,:,2,idx);
    xHat = HTestHat(:,:,1,idx) + 1i*HTestHat(:,:,2,idx);

    % Calculate correlation
    c1 = sqrt(sum(conj(x).*x,'all'));
    c2 = sqrt(sum(conj(xHat).*xHat,'all'));
    aa = abs(sum(conj(x).*xHat,'all'));
    correlation(idx) = aa/(c1*c2);

    % Calculate NMSE
    mseTmp = mean(abs(x-xHat).^2,'all');
    nmse(idx) = 10*log10(mseTmp/mean(abs(x).^2,'all'));

end

figure
tiledlayout(2,1)
nexttile
histogram(correlation,"Normalization","probability")
grid on
title(sprintf("Autoencoder Correlation (Mean Correlation = %1.5f)", ...
    mean(correlation)))
xlabel("Correlation"); ylabel("PDF")
nexttile
histogram(nmse,"Normalization","probability")
grid on
title(sprintf("Autoencoder NMSE (Mean NMSE = %1.2f dB)",mean(nmse)))
xlabel("NMSE (dB)"); ylabel("PDF")
sgtitle(sprintf('Tapped Delay Line %s Scenario',channel.DelayProfile));

%% Test network
[encNet,decNet] =  Autoencoder.SplitEncoderDecoder(trainedNet,'Enc_Sig');

% Generate new data
rng(1)
numFrames2 = 100;
numSamples2 = ceil(numFrames2*nRx*aenOptions.NumSlotsPerFrame);
[Hest,~,aenOptions] = Autoencoder.DataGeneration.GenerateData( ...
    numSamples2,simParams.Carrier,channel,aenOptions);

% Encode
aenOptions.Normalization = true;
codeword = Autoencoder.Encode(encNet,Hest,aenOptions);

% Decode
aenOptions.NumSubcarriers = simParams.Carrier.NSizeGrid*12;
Hhat = Autoencoder.Decode(decNet,codeword,aenOptions);

%% Calculate correlation and NMSE for end to end CSI feedback system and plot
Havg = squeeze(mean(Hest,2));
correlationE2E = zeros(nRx,numFrames2);
nmseE2E = correlationE2E;
for rxIdx = 1:nRx
    for frIdx = 1:numFrames2
        x = Havg(:,rxIdx,:,frIdx);
        xHat = Hhat(:,rxIdx,:,frIdx);

        % Compute correlation
        c1 = sqrt(sum(conj(x).*x,'all'));
        c2 = sqrt(sum(conj(xHat).*xHat,'all'));
        aa = abs(sum(conj(x).*xHat,'all'));
        correlationE2E(rxIdx,frIdx) = aa/(c1*c2);

        % Compute NMSE
        mseTmp = mean(abs(x-xHat).^2,'all');
        nmseE2E(rxIdx,frIdx) = 10*log10(mseTmp/mean(abs(x).^2,'all'));
    end
end
figure
tiledlayout(2,1)
nexttile
histogram(correlationE2E,"Normalization","probability")
grid on
title(sprintf("End-to-End Correlation (Mean Correlation = %1.5f)", ...
  mean(correlationE2E,'all')))
xlabel("Correlation"); ylabel("PDF")
nexttile
histogram(nmseE2E,"Normalization","probability")
grid on
title(sprintf("End-to-End NMSE (Mean NMSE = %1.2f dB)", ...
  mean(nmseE2E,'all')))
xlabel("NMSE (dB)"); ylabel("PDF")
sgtitle(sprintf('Tapped Delay Line %s Scenario',channel.DelayProfile));

%% Effects of Quantization
maxVal = 1;
minVal = -1;
idxBits = 1;
nBitsVec = 2:10;
correlationQ = zeros(nRx,numFrames2,length(nBitsVec));
nmseQ = zeros(nRx,numFrames2,length(nBitsVec));
for numBits = nBitsVec
    disp("Running for " + numBits + " bit quantization")

    % Quantize between 0:2^n-1 to get bits
    qCodeword = uencode(codeword*2 - 1, numBits);

    % Get back the floating point, quantized numbers
    codewordRx = (double(udecode(qCodeword,numBits)) + 1)/2;
    Hhat2 = Autoencoder.Decode(decNet,codewordRx,aenOptions);
    for rxIdx = 1:nRx
        for frIdx = 1:numFrames2
            x = Havg(:,rxIdx,:,frIdx);
            xHat = Hhat2(:,rxIdx,:,frIdx);

            % Compute correlation
            c1 = sqrt(sum(conj(x).*x,'all'));
            c2 = sqrt(sum(conj(xHat).*xHat,'all'));
            aa = abs(sum(conj(x).*xHat,'all'));
            correlationQ(rxIdx,frIdx,idxBits) = aa/(c1*c2);

            % Compute NMSE
            mseTmp = mean(abs(x-xHat).^2,'all');
            nmseQ(rxIdx,frIdx,idxBits) = 10*log10(mseTmp/mean(abs(x).^2,'all'));
        end
    end
    idxBits = idxBits + 1;
end

figure
tiledlayout(2,1)
nexttile
plot(nBitsVec,squeeze(mean(correlationQ,[1 2])),'*-')
title("Correlation (Codeword-" + size(codeword,3) + ")")
xlabel("Number of Quantization Bits"); ylabel("CorrelationQ")
grid on
nexttile
plot(nBitsVec,squeeze(mean(nmseQ,[1 2])),'*-')
title("NMSE (Codeword-" + size(codeword,3) + ")")
xlabel("Number of Quantization Bits"); ylabel("NMSE (dB)")
grid on
sgtitle(sprintf('Tapped Delay Line %s Scenario',channel.DelayProfile));
