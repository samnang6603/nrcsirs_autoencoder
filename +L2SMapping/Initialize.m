function L2SM = Initialize(carrier)
% Link-to-System (L2S) Mapping interface initializer
% Adapted from MATLAB L2SM internal functions
%
% L2S Mapping is the translator between PHY-level "link" simulations 
% (highly-detailed, bit-by-bit RF chain) and large-scale "system" 
% simulations (whole networks, coverage, mobility, interference). 
%
% In other words, it is an abstraction interface between core PHY layer, 
% aka Link-Level Simulator (LLS) to System-Level Simulator (SLS) network 
% transport traffics.
%
% Instead of running a full PHY sim every time a user clock in the network,
% we precompute link-level results (e.g., BLER vs SINR curves for each MCS)
% and then abstract them into a simpler model that system-level sims can 
% use.
%
% How is it used?
% ---------------
% 1. Compute Effective SINR: combine all the per-subcarrier/per-layer 
%    SINRs into one "effective SINR" value using:
%    - Effective Exponential SINR Mapping (EESM)
%    - Mutual Information Effective SINR Mapping (MIESM)
% 2. Map Effective SINR to BLER: look it up against precomputed curves from
%    link-level simulations (aka the "truth").
% 3. Use BLER to estimate throughput, latency, scheduling success, etc. 
%    at the system level.
%
% Why do we need it?
% ------------------
% - Link-level sims = accurate but slow (simulate hours of OFDM symbols).
% - System-level sims = scalable but too abstract (they don't "see" coding,
%   modulation, HARQ, or beamforming on their own).
% - L2S mapping is the middle ground: fast enough for system-level, 
%   detailed enough to stay faithful to real PHY behavior.
%

L2SM = struct();

% Carrier
L2SM.Carrier = carrier;
L2SM.GridSize = [];
L2SM.IndicesInfo = [];

% Received Bit Information Rate (RBIR)
L2SM.RBIR.Alpha = 1;
L2SM.RBIR.Beta  = 1;

% CQI Information
L2SM.CQI.Table = [];
L2SM.CQI.XOverhead = [];
L2SM.CQI.TableRowCombos = [];
L2SM.CQI.TableQms = [];
L2SM.CQI.TableModulations = [];
L2SM.CQI.Qm = [];
L2SM.CQI.G = [];
L2SM.CQI.TransportBlockSize = [];
L2SM.CQI.C = [];
L2SM.CQI.NBuffer = [];
L2SM.CQI.EffectiveCodeRate = [];
L2SM.CQI.DLSCHInfo = [];

% Limited Buffer Rate Matching (LBRM)
L2SM.Nref = [];

% Split code block
L2SM.SplitCodeBlock = false;

end