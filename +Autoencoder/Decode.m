function Hhat = Decode(dec,codeword,options)
%DecodeCSI decodes CSI using the encoder network 'enc' on codeword with 
%   options from auotencoder option 'opt'.

[nRx,N,codewordLen] = size(codeword);
nsc = options.NumSubcarriers;

tmp = zeros(1,codewordLen);
placeholder = predict(dec,tmp);
[~,nTx,~] = size(placeholder);

Hhat = zeros(nsc,nRx,nTx,N);
HTruncHat = zeros(options.MaxDelay,nTx,2,nRx*N,like=1i);

for rxIdx = 1:nRx
    for nIdx = 1:N
        % Decode
        tmp = codeword(rxIdx,nIdx,:);
        tmpInfer = predict(dec,tmp(:).');
        HTruncHat(:,:,:,(nIdx-1)*nRx + rxIdx) = tmpInfer;
        Hhat(:,rxIdx,:,nIdx) = Autoencoder.DataGeneration.PostprocessChannelResponse(tmpInfer,options);
    end
end

end