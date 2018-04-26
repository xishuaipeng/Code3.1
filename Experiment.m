%%%%%%%%%%%%%%%%%%%%%Experiment 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc;
close all;
addpath './OBDExtraction';
addpath './DataClass';
addpath './Feature'
addpath './Feature/VGG16'
addpath './Net'
videoName = {'112_07172017','118_07182017', '106_07142017','023','028'};
data = [];
label = [];
tic;
%     obj.localParameters.addParamValue('start_padding', 0.1);
%     obj.localParameters.addParamValue('last_padding', 0.02);
%       obj.localParameters.addParamValue('gene_feature', {'Curvature','VGGObject','Heading','Speed'});
for i =1:length(videoName)
    sprintf('%s is processing!',videoName{i})
    Mdata = Maneuverdata('data_id', videoName{i},...
        'input_path','D:\TRI\NXR\input',...
        'out_dir','D:\TRI\XR3.0\Pro',...
        'start_padding',0.1,...
         'last_padding',0.1,...
         'gene_feature', {'Curvature','VGGScene','VGGObject','Heading','Speed'},...
         'load_feature', {'Curvature','VGGScene','Heading','Speed'}...
        );
    [seqX,seqY] = Mdata.EventSequence();
    data = [data;seqX];
    label = [label;seqY];
end
toc;
%%%%%%%%
[trainD,trainL,testD, testL,testIndex] = selectData(data,label.eventType,0.7);
net = ManeuversNet(trainD, trainL,100,1000);
y = classify(net,testD );
[testAna, accuracy] = generateAccuracy(double(testL)-1,double(y)-1)
corIndex= find( (double(testL) - double(y)) == 0 );
corInf = label(corIndex,:);
fprintf('sequence end to event end(time): %f \n', mean(corInf.seqEnd2eventEnd_time) );
fprintf('sequence end to event end(distance): %f \n', mean(corInf.seqEnd2eventEnd_distance) );
fprintf('sequence end to begin end(time): %f \n', mean(corInf.seqEnd2eventBeg_time) );
fprintf('sequence end to begin end(distance): %f \n', mean(corInf.seqEnd2eventBeg_distance) );






