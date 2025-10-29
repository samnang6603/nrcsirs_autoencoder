function csiFeedbackOpts = ConfigureCSIFeedbackOptions(simParams,snrIdx)

csiFeedbackOpts.CSIReportConfig = simParams.CSIReportConfig;
csiFeedbackOpts.DMRSConfig = simParams.PDSCH.DMRS;
csiFeedbackOpts.PerfectChannelEstimator = simParams.PerfectChannelEstimator;

% if strcmpi(simParameters.CSIReportMode,"AI CSI compression")
%     % Copy additional link adaptation configuration for AI CSI compression mode
%     csiFeedbackOpts.AINetworkFilename = simParameters.AINetworkFilename;
% 
%     % Download and extract a pretrained CSI network for AI CSI compression mode
%     displayProgress = (snrIdx==1);
%     helperCSINetDownloadData(displayProgress);
% end

end