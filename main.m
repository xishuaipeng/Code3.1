clear,clc;
close all;
addpath './OBDExtraction';
addpath './DataClass';
addpath './Feature'
addpath './Feature/VGG16'
addpath './Net'%106_07142017%'023','106_07142017', '028', '112_07172017','118_07182017','ID002_T001','ID002_T002','ID002_T003','ID002_T004','ID002_T010','ID002_T015'
videoName = {  '112_07172017','118_07182017', '106_07142017','023','028'};
% '106_07142017','023','028', ...
%     'ID001_T001','ID001_T002','ID001_T003','ID001_T004','ID001_T006','ID001_T007', ...
%     'ID001_T008','ID001_T009','ID001_T010','ID001_T011','ID001_T012','ID001_T013',...
%     'ID001_T014','ID001_T015', 'ID001_T017',...
%     'ID002_T001','ID002_T002','ID002_T003','ID002_T004','ID002_T010', ...
%     'ID002_T012','ID002_T013','ID002_T015','ID002_T016','ID002_T017', ...
%     'ID002_T018','ID002_T019','ID002_T020','ID002_T021', ...
%     'ID003_T001','ID003_T002','ID003_T003', ...
%     'ID003_T004','ID003_T010','ID003_T014'
%'106_07142017', '112_07172017','118_07182017','023','028','ID002_T001','ID002_T002','ID002_T003','ID002_T004','ID002_T010','ID002_T015','ID002_T016','ID002_T017','ID002_T018','ID002_T019','ID002_T020','ID002_T021'};%'ID002_T001','ID002_T002','ID002_T003',}%'106_07142017, '023','028','112_07172017','118_07182017','ID002_T001','ID002_T002','ID002_T003','ID002_T004','ID002_T010','ID002_T015','ID002_T016','ID002_T017','ID002_T018','ID002_T019','ID002_T020','ID002_T021'};%'ID002_T001','ID002_T002','ID002_T003',
% 'ID001_T016', 'ID002_T009',bad trip
% videoName = {'ID001_T001','ID001_T010','ID001_T012','ID001_T013','ID001_T015', ...
%    'ID001_T017','ID003_T001','ID003_T002','ID003_T003', ...
%     'ID003_T004','ID003_T010','ID003_T014'};
data = [];
label = [];
minX = 1e10;
maxX = 0;
sumX = 0;
squareX = 0;
numX = 0;
tic;
for i =1:length(videoName)
    sprintf('%s is processing!',videoName{i})
    Mdata = Maneuverdata('data_id', videoName{i},'input_path','D:\TRI\NXR\input','out_dir','D:\TRI\XR3.0\Pro');
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
% netList = dir('./model');
% for i = 3:100:length(netList)
%     netTest = ['./model/',netList(i).name];
%     load(netTest)
%     y = classify(net,testD );
%     [testAna, accuracy] = generateAccuracy(double(testL)-1,double(y)-1);
%     sprintf('%s:%f',netList(i).name,accuracy),
% end
%accuracy =  sum(y == testL)/numel(testL)
%Extra Test
% load('process.mat');
% Mdata = Maneuverdata('ID002_T021');
% Mdata = Mdata.GenerateSequence();
% X = Mdata.X;
% X = cellfun(@(x) (x-EX)./(varX), X,'UniformOutput',false);
% y = classify(net,X );
% [testAna, accuracy] = generateAccuracy(Mdata.Y,double(y)-1)
% save('net.mat','net');





