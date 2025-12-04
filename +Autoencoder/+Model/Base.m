function aen = Base
%CreateAutoencoder creates autoencoder network based on a prototype from
%   MATLAB example.
%
%
% Reference:
% 1. https://www.mathworks.com/help/releases/R2025a/5g/ug/csi-feedback-with-autoencoders.html

encoderLayers = [
    imageInputLayer([28, 2, 2],Normalization='none',Name='Enc_Input')
    convolution2dLayer([3, 3],2,Padding='same',Name='Enc_Conv')
    batchNormalizationLayer(Epsilon=0.001,Name='Enc_BN')
    leakyReluLayer(0.3,Name='Enc_LeakyRelu')
    flattenLayer(Name='Enc_Flatten')
    fullyConnectedLayer(64,Name='Enc_FCN')
    sigmoidLayer(Name='Enc_Sig')
    ];

%aen = addLayers(aen,encoderLayers);

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
end