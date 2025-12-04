function Hestpp = PreprocessChannelResponse(Hest,options)
%preprocessChannelResponse Prepares the channel estimate 'Hest' by
%   decimating via 2D-FFT and 2D-IFFT. The preprocessed channel estimate
%   has maxDelay samples and is a real-valued tensor.
%   Hestpp is maxDelay x nTx x 2(I/Q) x nRx x Nslots, where 
%   Nslots = numSyms/14
%
%   This procedure is similar to Principle Components Analysis (PCA)
%

maxDelay = options.MaxDelay;
[numSubcarrier,numSyms,nRx,nTx] = size(Hest);
% Compute center of delay spread band. Keep only the taps that matter and
% discard channel where magnitude is close to 0. This is delay-domain
% sparsity exploitation.
midpt = floor(numSubcarrier/2);
lb = midpt - (numSubcarrier - maxDelay)/2 + 1; % lower bound
ub = midpt + (numSubcarrier - maxDelay)/2; % upper bound

% Because fft/fft2 in MATLAB arranges in natural frequency order (not
% 0 Hz-centered fftshifted, the lower bound represents "positive frequency"
% and the upper bound represents "negative frequency".

%{
The delay-domain truncation uses bounds lb and ub computed around midpt
(the FFT center of the subcarrier axis). This may appear counterintuitive,
because physical multipath delays are non-negative. However, when we obtain
delay taps by taking an FFT across subcarriers, the discrete Fourier 
transform (DFT) represents delays as complex exponentials in the frequency 
index k. A delay l introduces a linear phase ramp exp(-j2*pi*k*l/N). The 
DFT of this ramp produces a spectral line whose index is interpreted as a 
"frequency bin".

Due to the periodic nature of the DFT, delay taps appear in BOTH the
positive-frequency and negative-frequency halves of the FFT output:

  - Small delays -> low-frequency bins near index 0
  - Large delays -> bins near index N (wrapped-around negative frequencies)

Thus, the effective delay spread (maxDelay taps) manifests as TWO clusters:
one at the beginning of the FFT vector and one at the end. The region 
between these two clusters contains only noise and should be discarded.

The indices lb and ub mark this central noise-only region to be removed. The
kept taps reside in:
    1 : (lb-1)       and       (ub+1) : N

which together contain exactly maxDelay useful delay samples.

FFT delay axis (no fftshift):

FFT index-> 0                          midpt                            N-1
            |------------------------|discard|----------------------------|
            |  positive frequencies  |region |    negative frequencies    |
            |   (small delays)       |.......|    (wrapped large delays)  |
            |<----- useful taps ---->|.......|<------ useful taps ------->|
             ^ left keep region               ^ right keep region
                    (lb-1)                         (ub+1)

The central block between lb and ub contains negligible channel energy and 
is therefore truncated for compression. During reconstruction, zeros are 
inserted back into this region to restore the full FFT size before applying
the IFFT.

EXAMPLE:

For N = 624 subcarriers:
    midpt = floor(N/2) = 312
    maxDelay = 28
    lb = midpt - (N - maxDelay)/2 + 1 = 15
    ub = midpt + (N - maxDelay)/2     = 610

This means:
    - We KEEP bins: [ 1 : 14 ] and [ 611 : 624 ]  --> total = 28 taps
    - We DISCARD bins: [ 15 : 610 ]               --> noise-only region

Why this happens:
-----------------
An FFT across subcarriers maps multipath delays into *DFT frequency bins*.
Because the DFT is cyclic, delay taps appear split between the beginning
("positive frequencies") and end ("negative frequencies") of the FFT 
output. Both correspond to physical delays; the division is purely an 
artifact of DFT indexing (no fftshift used).

ASCII diagram (to scale in concept, not exact spacing):

0                               midpt (=312)                     N-1 (=623)
|----------------------------| discard region |---------------------------|
|     positive frequencies   |    (noise)     |     negative frequencies  |
|       (small delays)       |................|    (wrapped large delays) |
|<----- useful taps -------->|................|<------- useful taps ----->|
keep: 1..14                      ^ discard ^            keep: 611..624
                                      |
                               lb = 15, ub = 610

Thus, the useful delay spread (maxDelay = 28 taps) manifests as:
    - 14 taps at the beginning of the FFT output
    - 14 taps at the end (wrapped-around)
These together represent the true CIR support.

During preprocessing:
    - We discard the middle region [lb : ub] that contains no meaningful 
      energy (decimation).

During postprocessing:
    - We reinsert zeros at [lb : ub] to reconstruct the original FFT length
      before applying the 2D IFFT (interpolation).

This behavior is equivalent to truncating the channel's delay spread around
0 delay, but expressed in the DFT's natural indexing order.

%}

% Average over symbols
symsPerSlot = 14;
Nslots = numSyms/symsPerSlot;
HTmp = reshape(Hest,numSubcarrier,Nslots,symsPerSlot,nRx,nTx);
% Remove fast fading noise and keeps long-term multipath structure
H = mean(HTmp,3); % 3rd dim is symsPerSlot dim
H = permute(H, [1, 5, 4, 2, 3]); % rearrange dims

% Pre-allocate
Hestpp = zeros(maxDelay,nTx,2,nRx,Nslots,like=1i);

% Per Rx antenna, decimate over subcarriers using 2D-FFT
for slotIdx = 1:Nslots
    % 2D-FFT: go from freq/spatial(antenna) to delay/angle(spatial dir) 
    % domain
    Hfft = fft2(H(:,:,:,slotIdx));

    % Truncate (explained at lb, ub calculation section)
    HTruncTmp = Hfft([1:lb-1, ub+1:end],:,:);

    switch options.DataDomain
        case 'Frequency-Spatial'
            % Transform back if freq/spatial is desired
            HTrunc = ifft2(HTruncTmp);
        otherwise
            % Keep as is iff delay/spatial is desired
            HTrunc = HTruncTmp;
    end

    % Extract Real (I) and Imag (Q) components
    HI = real(HTrunc); % I
    HQ = imag(HTrunc); % Q

    if options.Normalization
        avgVal = options.AverageValue;
        stdVal = options.StandardDeviationValue;
        targetStd = options.TargetStandardDeviationValue;
        % Normalize using standard Z-score and rescaling std
        Hestpp(:,:,1,:,slotIdx) = ((HI - avgVal)/stdVal)*targetStd + 0.5;
        Hestpp(:,:,2,:,slotIdx) = ((HQ - avgVal)/stdVal)*targetStd + 0.5;
    else
        % Keep as is
        Hestpp(:,:,1,:,slotIdx) = HI;
        Hestpp(:,:,2,:,slotIdx) = HQ;
    end
end
