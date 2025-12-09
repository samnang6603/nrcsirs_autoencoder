function [net,autoEncOpt,encEndLayer] = LoadPretrainedNetwork(csiReportConfig)
%loadPretrainedAutoencoder Load pretrained autoencoder from mathworks
%   example. The example can be found here:
%   https://www.mathworks.com/help/releases/R2025a/5g/ug/nr-pdsch-throughput-using-csi-feedback.html

modelName = csiReportConfig.Autoencoder.ModelName;
if ispc
    modelFilePath = '+Autoencoder\Pretrained_Models\';
else
    modelFilePath = '+Autoencoder/Pretrained_Models/';
end
if isfolder(modelFilePath)
    load([modelFilePath,modelName],'net','autoEncOpt')
else
    error('Cannot locate Pretrained_Models directory')
end
% Add some parameters to pretrained network from Mathworks to fit custom
% model
autoEncOpt.DataDomain = 'Frequency-Spatial';
autoEncOpt.AverageValue = autoEncOpt.MeanVal;
autoEncOpt.StandardDeviationValue = autoEncOpt.StdValue;
autoEncOpt.TargetStandardDeviationValue = autoEncOpt.TargetSTDValue;

% Known encoder end layer
encEndLayer = "Enc_Sigmoid";
end