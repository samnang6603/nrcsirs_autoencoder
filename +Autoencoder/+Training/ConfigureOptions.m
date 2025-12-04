function [trainOpts,lossFcn] = ConfigureOptions(solverName,valData,miniBatchSize,maxEpochs)
%CreateTraining Outputs a training options for general deep learning module

if nargin < 4
    maxEpochs = 1000;
end

trainOpts = trainingOptions(solverName, ...
    InitialLearnRate=0.01, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropPeriod=156, ...
    LearnRateDropFactor=0.5916, ...
    Epsilon=1e-7, ...
    MaxEpochs=maxEpochs, ...
    MiniBatchSize=miniBatchSize, ...
    Shuffle="every-epoch", ...
    ValidationData=valData, ...
    ValidationFrequency=20, ...
    Metrics="rmse", ...
    Verbose=true, ...
    ValidationPatience=20, ...
    OutputNetwork="best-validation-loss", ...
    ExecutionEnvironment="auto", ...
    Plots='training-progress');

lossFcn = @(x,t) nmseLossdB(x,t);

end

%% Local Helper Fcn
function loss = nmseLossdB(x,xHat)
%nmseLossdB NMSE loss in dB
in = complex(x(:,:,1,:),x(:,:,2,:));
out = complex(xHat(:,:,1,:),xHat(:,:,2,:));
nmsePerObservation = helperNMSE(in,out);
loss = mean(nmsePerObservation);
end

function nmse = helperNMSE(in,out,varargin)
%helperNMSE Normalized mean square error
%   NMSE = helperNMSE(X,Y) returns the normalized mean square
%   error (NMSE) between X and Y. X and Y must be complex valued column
%   vectors or 2-by-N array where the real and imaginary parts of the
%   signal is carried on the first and second rows, respectively.
%
%   NMSE = helperNMSE(X,Y,FORMAT) returns the NMSE between X and Y in the
%   specified format, which can be "dB" or "linear". The default is
%   "dB". 

%   Copyright 2023-2024 The MathWorks, Inc.

if ~isreal(in)
  mse = mean(abs(in-out).^2,'all');
  nmse = mse/mean(abs(in).^2,'all');
else
  mse = mean(sum((in-out).^2,1),'all');
  nmse = mse/mean(sum(in.^2,1),'all');
end

if nargin > 2
  format = varargin{1};
else
  format = "dB";
end
if strcmp(format,"dB")
  nmse = 10*log10(nmse);
end
end