% End-to-end CSI Feedback using AEN

[fname, fpath, ~] = uigetfile("*.mat",MultiSelect="on");
load([fpath,fname],'trainedNet');

[encNet,decNet] =  Autoencoder.SplitEncoderDecoder(trainedNet,'Enc_Sig');

figure,
subplot(131)
plot(trainedNet), title('Autoencoder')
subplot(132)
plot(encNet), title('Encoder')
subplot(133)
plot(decNet), title('Decoder')

