% Prototype of Autoencoder for CSI-RS Report Feedback
% Reference:
% 1. https://www.mathworks.com/help/releases/R2025a/5g/ug/csi-feedback-with-autoencoders.html
% 2. Zimaglia, Elisa, Daniel G. Riviello, Roberto Garello, and Roberto Fantini.
%   "A Novel Deep Learning Approach to CSI Feedback Reporting for NR 5G Cellular Systems." 
%   In 2020 IEEE Microwave Theory and Techniques in Wireless Communications (MTTW), 
%   47–52. Riga, Latvia: IEEE, 2020.
%   https://doi.org/10.1109/MTTW51045.2020.9245055.

clear, clc

trainingToggle = true;

aen = dlnetwork();

encoderLayers = [
    imageInputLayer([7, 2, 2],Normalization='none',Name='Enc_Input')
    convolution2dLayer([3, 3],2,Padding='same',Name='Enc_Conv')
    batchNormalizationLayer(Epsilon=0.001,Name='Enc_BN')
    leakyReluLayer(0.3,Name='Enc_LeakyRelu')
    flattenLayer(Name='Enc_Flatten')
    fullyConnectedLayer(64,Name='Enc_FCN')
    sigmoidLayer(Name='Enc_Sig')
    ];

aen = addLayers(aen,encoderLayers);

decoderStage1Layers = [
    fullyConnectedLayer(28,Name='Dec_FCN')
    functionLayer(@(x)dlarray(reshape(x,7,2,2,[]),'SSCB'),Formattable=true,...
                  Acceleratable=true,Name='Dec_Reshape_S1')

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
aen = connectLayers(aen,'Enc_Sig','Dec_FCN');
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

if trainingToggle
    % Carrier
    simParams.Carrier = nrCarrierConfig;

    % Antenna Configurations
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

    % Channel Delay Profile
    simParams.DelayProfile = 'CDL-C';   % 'CDL-' or 'TDL-'
    simParams.DelaySpread = 300e-9;     % s
    simParams.MaximumDopplerShift = 5;  % Hz

    % Create Channel
    rng(0);
    % channel = Channel.CreateChannel(simParams);
    % channel.ChannelFiltering = false;
    % channel.ChannelResponseOutput = 'path-gains';
    % channelInfo = info(channel);
    channel = nrCDLChannel;
    channel.DelayProfile = 'CDL-C';
    channel.DelaySpread = 300e-9;       % s
    channel.MaximumDopplerShift = 5;   % Hz
    channel.RandomStream = "Global stream";
    channel.TransmitAntennaArray.Size = [1 1 2 1 1];
    channel.ReceiveAntennaArray.Size = [1 1 2 1 1];
    channel.ChannelFiltering = false;           % No filtering

    % Autoencoder Options
    aenOptions.DataDomain = 'Frequency-Spatial';
    aenOptions.TruncationFactor = 10;
    aenOptions.ZeroTimingOffset = true;
    aenOptions.SaveData = false;
    aenOptions.Verbose = true;

    % Generate channel response data
    numSamples = 1500;
    [~,HTargetpp,aenOptions] = Autoencoder.DataGeneration.GenerateData(numSamples,...
        simParams.Carrier,channel,aenOptions);

    % Each frame has a data for nRx 
    [maxDelay,nTx,NIQ,nRx,Nframes] = size(HTargetpp);

    % Combine frames and nRx
    HTargetpp = reshape(HTargetpp,maxDelay,nTx,NIQ,nRx*Nframes);

    % Calculate the average and standard deviation to normalize the data
    avgVal = mean(HTargetpp,'all');
    stdVal = std(HTargetpp,[],'all');

    % Separate data into training, validation and testing
    numDataSet = size(HTargetpp,4);
    numTrainDS = 0.6*numDataSet;
    numValidDS = 0.15*numDataSet;
    numTestDS  = 0.25*numDataSet;

    % Target Std of value of 0.0212, which will restrict data to range
    % [-0.5 , 0.5]
    trainInd = 1:numTrainDS;
    validInd = numTrainDS + (1:numValidDS);
    testInd  = numTrainDS + numValidDS + (1:numTestDS);
    HTrain = (HTargetpp(:,:,:,trainInd) - avgVal)/stdVal*0.0212 + 0.5;
    HValid = (HTargetpp(:,:,:,validInd) - avgVal)/stdVal*0.0212 + 0.5;
    HTest  = (HTargetpp(:,:,:,testInd)  - avgVal)/stdVal*0.0212 + 0.5;

    % Update aen options
    aenOptions.AvgVal = avgVal;
    aenOptions.StdVal = stdVal;
    aenOptions.TargetStd = 0.0212;

    
    




    




end


