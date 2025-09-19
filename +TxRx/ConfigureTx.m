function [carrier,encDLSCH,pdsch,pdschextra,csirs,wtx] = ConfigureTx(simParams)
% Get channel and signal-level params create DL encoders and precoding
% matrix

carrier = simParams.Carriers;
pdsch   = simParams.PDSCH;
pdschextra = simParams.PDSCHExtension;
csirs   = simParams.CSIRS;

% DLSCH Encoder
encDLSCH = nrDLSCH;

% Precoding matrix initialization
wtx = 1;

% Determine xOverhead: maximum allowed CSI-RS RE overhead for the given BWP,
% rounded to the nearest valid bucket. This ensures the CSI-RS allocation 
% stays within 3GPP spec limits and does not affect PDSCH decoding or other
% critical channels. 
% The same concept applies to SRS and PUCCH/PUSCH allocations.
if iempty(pdschextra.xOverhead)
    [~,csirsInfo] = nrCSIRSIndices(carrier,csirs);

    % Get Frequency-domain and time-domain locations of the lowest resource 
    % elements corresponding to all code division multiplexing (CDM) groups
    % Basically, the "starting points" of each CDM group
    KBarLBarLen = length(csirsInfo.KBarLBar{1}); 

    % CDM Group: A set of resource elements (REs) that are multiplexed via 
    % code division (different sequences / polarizations).
    % Each CDM group has a "base" location KBarLBar and then small offsets
    % KPrime, LPrime to locate each RE inside that group

    % Get Frequency-domain indexing (offset) within a CDM group
    KPrimeLen = length(csirsInfo.KPrime{1});

    % Get Time-domain indexing (offset) within a CDM group
    LPrimeLen = length(csirsInfo.LPrime{1});

    % The total REs for the CSI-RS
    csirsRE = KBarLBarLen*KPrimeLen*LPrimeLen;

    % TS 38.214 5.1.3.2 Page 29
    xOverheadBuckets  = [0, 6, 12, 18]; % From xOverhead
    partition = [0, 6, 12]; % range endpoints

    % Calculate xOverhead as a quantized value from the xOverheadBucket
    [~,xOverhead] = quantiz(csirsRE,partition,xOverheadBuckets);
    pdschextra.xOverhead = xOverhead;

    if csirsRE > xOverhead
        warning("CSI-RS RE usage is hgiher than the maximum allowed of 18..." + ...
            ". A decoding error might occur at PDSCH/PUSCH")
    end


end

end