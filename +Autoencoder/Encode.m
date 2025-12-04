function codeword = Encode(enc,Hest,options)
%EncodeCSI encodes CSI using the encoder network 'enc' on Hest with options
% from auotencoder option 'opt'.

[~,~,nRx,nTx,D] = size(Hest);

tmp = zeros(options.MaxDelay,nTx,2);

placeholder = predict(enc,tmp);

codeword = zeros(nRx,D,length(placeholder));

for idx = 1:D
    Hpp = Autoencoder.DataGeneration.PreprocessChannelResponse(Hest(:,:,:,:,idx),options);
    tmpInfer = predict(enc,Hpp);
    codeword(:,idx,:) = tmpInfer;
end


end