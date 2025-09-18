function reportConfig = ConfigureCodebookParameters(reportConfig)
%ConfigureCodebookParameters Configure panel for type I single-panel
%codebook. 
% 
%
% Notes:
%   - NumCSIRSPorts (from CSI-RS config) must equal 2*N1*N2.
%   - nHEstCSIRSPorts represents the number of CSI-RS ports in the 
%     estimated channel matrix H (K×L×nRxAnts×nHEstCSIRSPorts).
%

% Supported panel-configuration table from TS 38.214 Table 5.2.2.2.2.1-2
% N1: The number of antenna elements in the horizontal direction
% N2: The number of antenna elements in the vertical direction
% O1: The DFT oversampling factor for antenna elements in the horizontal direction
% O2: The DFT oversampling factor for antenna elements in the vertical direction
% type1SinglePanelConfig = [...
%     % Number of CSI-RS Port   (N1,N2)  (O1,O2)
%                           4  ,  2,1   ,  4,1
%                           8  ,  2,2   ,  4,4
%                           8  ,  4,1   ,  4,1
%                          12  ,  3,2   ,  4,4
%                          12  ,  6,1   ,  4,1
%                          16  ,  4,2   ,  4,4
%                          16  ,  8,1   ,  4,1
%                          24  ,  4,3   ,  4,4
%                          24  ,  6,2   ,  4,4
%                          24  , 12,1   ,  4,1
%                          32  ,  4,4   ,  4,4
%                          32  ,  8,2   ,  4,4
%                          32  , 16,1   ,  4,1
%                          ]
portcol = [4,8,8,12,12,16,16,24,24,24,32,32,32].';
N1col = [2,2,4,3,6,4,8,4,6,12,4,8,16].';
N2col = [1,2,1,2,1,2,1,3,2,1,4,2,1].';
O1col = 4*ones(13,1);
O2col = [1,4,1,4,1,4,1,4,4,1,4,4,1].';
Type1SinglePanelConfigTable = table(portcol,N1col,N2col,O1col,O2col,...
    VariableNames=["Port","N1","N2","O1","O2"]);

% Get N1 and N2 from Panel Dimensions
N1 = reportConfig.PanelDimensions(1); % horizontal
N2 = reportConfig.PanelDimensions(2); % vertical

% PCSIRS (CSI-RS port)
Pcsirs = N1*N2*2;
reportConfig.Pcsirs = Pcsirs;

% Find the configuration index corresponding to N1 and N2
N1ConfigIdx = Type1SinglePanelConfigTable.N1 == N1;
N2ConfigIdx = Type1SinglePanelConfigTable.N2 == N2;
configIdx = find(N1ConfigIdx & N2ConfigIdx,1); % return the first index that matches both

% Then get the oversampling factors
O1 = Type1SinglePanelConfigTable.O1(configIdx);
O2 = Type1SinglePanelConfigTable.O2(configIdx);

% Update reportConfig field
reportConfig.OverSamplingFactors = [O1,O2];

% Codebook's length is the product of N1,N2,O1,O2 (total number of precoder
% entries) and the codebookSubsetRestriction is a vector of ones whose 
% length is the codebook's length. It is a placeholder representing
% default/unconstrained case.
codebookLen = N1*N2*O1*O2;
reportConfig.CodebookSubsetRestriction = ones(1,codebookLen);

% i2Restriction: i2 refinement index restriction within the (i1,i2) chosen
% precoder. i2 is the bit sequence defined on page 66 of TS 38.214
reportConfig.i2Restriction = ones(1,16);

end