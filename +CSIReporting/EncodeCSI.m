function csiReport = EncodeCSI(carrier,csirs,Hest,nVar,csiFeedbackOpts)
%EncodeCSI Prepares and encodes CSI

csiReportConfig = csiFeedbackOpts.CSIReportConfig;
switch csiReportConfig.Mode
    case 'AUTOENCODER'
        % Using the encoder portion of the autoencoder to encode CSI
        csiReport = encodeViaAutoencoder(carrier,Hest,nVar,csiFeedbackOpts);
    case 'RI-PMI-CQI'
        % Using 3GPP TS 38.211/214 RI-PMI-CQI selection
        csiReport = encodeViaRIPMICQI(carrier,csirs,Hest,nVar,csiFeedbackOpts);
    otherwise
        error('Invalid CSI Report Mode')
end

end

%% Local Helper Fcns
function csiReport = encodeViaRIPMICQI(carrier,csirs,H,nVar,csiFeedbackOpts)
%encodeRIPMICQI Subroutine to encode RI-PMI-CQI reporting
% Adjust noise if practical channel
if ~csiFeedbackOpts.PerfectChannelEstimator
    % Empirical noise estimation error scaling from practical channel
    % estimator. 
    tmpdB = 2.0466; % Likely empirical according to MATLAB
    nVar = nVar*db2pow(tmpdB);
end

% Simulate the PMI-CQI for all possible rank and select the best rank
rankSearchCriterion = 'MaxSE'; % Force max spectral efficiency search
[ri,~,CQIPMICompParams] = CSIReporting.SelectRI(carrier,csirs,csiFeedbackOpts,H,nVar,rankSearchCriterion);

% If no available/possible rank, use rank 1
if isnan(ri)
    ri = 1;
end

% Use the chosen RI to select CQI
CQIPMICompParams.ThisRank = ri;
[cqi,pmi,~,pmiInfo] = CSIReporting.SelectCQI(CQIPMICompParams);

% Aggregate CSI results
csiReport = struct();
csiReport.CQI = cqi;
csiReport.RI  = ri;
csiReport.PMI = pmi;
csiReport.Precoder = pmiInfo.W;
csiReport.NSlot = CQIPMICompParams.Carrier.NSlot;

end

function csiReport = encodeViaAutoencoder(carrier,Hest,nVar,csiFeedbackOpts)
%selectCSIAutoencoder Encode CSI using autoencoder

% Load the encoder body
encNet  = csiFeedbackOpts.CSIReportConfig.Autoencoder.Encoder;
aenOpts = csiFeedbackOpts.CSIReportConfig.Autoencoder.Options;

% Using the encoder body, do inference on the estimated channel response
Hc = Autoencoder.Encode(encNet,Hest,aenOpts);

% Output structure
csiReport.H = Hc;
csiReport.nVar = nVar;
csiReport.NSlot = carrier.NSlot;
end