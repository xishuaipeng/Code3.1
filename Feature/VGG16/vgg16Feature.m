function [feature] = vgg16Feature(image,net)
curPath = mfilename('fullpath');
file  = strsplit( curPath,["/","\"]);
feature =0;
prototxtPath = fullfile(file{1:end-1},'deploy_vgg16_places365.prototxt');
caffemodelPath = fullfile(file{1:end-1},'vgg16_places365.caffemodel');
if isempty(net)
    net = importCaffeNetwork(prototxtPath,caffemodelPath);
end
sz = net.Layers(1).InputSize;
image = imresize(image,[sz(1),sz(2)]);
feature = activations(net,image,36);
feature = double(feature(:));
end

