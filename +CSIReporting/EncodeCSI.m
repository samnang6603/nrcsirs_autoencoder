function csiReport = EncodeCSI(carrier,csirs,Hest,nVar,csiFeedbackOpts)
%EncodeCSI Prepares and encodes CSI

switch csiFeedbackOpts.CSIReportMode
    case 'AI CSI compression'

    case 'RI-PMI-CQI'
        % Using 3GPP TS 38.211/214 RI-PMI-CQI Selection
        csiReport = selectCSI(carrier,csirs,H,nVar,csiFeedbackOpts);
    otherwise

end

end

function csiReport = selectCSI(carrier,csirs,H,nVar,csiFeedbackOpts)
%selectCSI Subroutine to run RI-PMI-CQI reporting
% Adjust noise if practical channel
if ~csiFeedbackOpts.PerfectChannelEstimator
    % Empirical noise estimation error scaling from practical channel
    % estimator. 
    tmpdB = 2.0466; % Likely empirical according to MATLAB
    nVar = nVar*db2pow(tmpdB);
end

% Simulate the PMI-CQI for all possible rank and select the best rank
rankSearchCriterion = 'MaxSE'; % Force max spectral efficiency search
ri = SelectRI(carrier,csirs,dmrsConfig,reportConfig,Hest,rankSearchCriterion);

end