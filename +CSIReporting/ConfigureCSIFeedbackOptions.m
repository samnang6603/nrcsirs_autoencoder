function csiFeedbackOpts = ConfigureCSIFeedbackOptions(simParams,snrIdx)

csiFeedbackOpts = struct();
csiFeedbackOpts.CSIReportMode = simParams.CSIReportMode;

csiFeedbackOpts.CSIReportPeriod = simParameters.CSIReportConfig.Period;
csiFeedbackOpts.CSIReportConfig = simParameters.CSIReportConfig;
csiFeedbackOpts.PerfectChannelEstimator = simParameters.PerfectChannelEstimator;
csiFeedbackOpts.DMRSConfig = simParameters.PDSCH.DMRS;

% if strcmpi(simParameters.CSIReportMode,"AI CSI compression")
%     % Copy additional link adaptation configuration for AI CSI compression mode
%     csiFeedbackOpts.AINetworkFilename = simParameters.AINetworkFilename;
% 
%     % Download and extract a pretrained CSI network for AI CSI compression mode
%     displayProgress = (snrIdx==1);
%     helperCSINetDownloadData(displayProgress);
% end

end