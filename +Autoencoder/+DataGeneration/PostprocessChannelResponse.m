function Hhat = PostprocessChannelResponse(Hpp,options)
%PostprocessChannelResponse applies post-processing for Hest by
%   interpolating via 2D-FFT and 2D-IFFT (reversal of preprocessing)

maxDelay = options.MaxDelay;
avgVal = options.AverageValue;
stdVal = options.StandardDeviationValue;
targetStdVal = options.TargetStandardDeviationValue;
nsc = options.NumSubcarriers;
nTx = size(Hpp,2);

midpt = floor(nsc/2); % midpoint
lb = midpt - (nsc - maxDelay)/2 + 1; % lower bound

HTmp = stdVal*((Hpp - 0.5)/targetStdVal) + avgVal;

% Combine real and imag from 3rd dim
HTmpCmplx = HTmp(:,:,1) + 1i*HTmp(:,:,2);
HTmpCmplxfft = fft2(HTmpCmplx);
HTmp2 = [HTmpCmplxfft(1:lb-1,:); % Insert up to the lower bound
         zeros((nsc - maxDelay),nTx); % pad in zeros (compensate for truncation in preprocessing)
         HTmpCmplxfft(lb:end,:)]; % Insert from the edge of lower bound to end
Hhat = ifft2(HTmp2);

end