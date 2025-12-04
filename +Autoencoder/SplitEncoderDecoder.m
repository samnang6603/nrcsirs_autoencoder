function [enc,dec] = SplitEncoderDecoder(aen,encEnd)
%SplitEncoderDecoder splits the autoencoder 'aen' into two separate encoder
%   and decoder parts

encEnd = string(encEnd);

encEndIdx = find({aen.Layers.Name} == encEnd);
decStart = aen.Layers(encEndIdx+1).Name;
decStartInputSize = aen.Layers(encEndIdx + 1).InputSize;

% Remove encoder layers from aen
allDecLayers = {aen.Layers(encEndIdx+1:end).Name};
enc = removeLayers(aen,allDecLayers);

% Remove decoder layers from aen
allEncLayers = {aen.Layers(1:encEndIdx).Name};
dec = removeLayers(aen,allEncLayers);

% Since decoder no longer has feature input layer, we will need to manually
% add to it
fL = featureInputLayer(decStartInputSize,Name="Dec_Input");
dec = addLayers(dec,fL);
dec = connectLayers(dec,"Dec_Input",decStart);

% Initialize
enc = enc.initialize();
dec = dec.initialize();

end