function [trainD,trainL,testD, testL] = SplitData(data,label,radio)
%SPLITDATA Summary of this function goes here
%   Detailed explanation goes here
minL = min(label);
maxL = max(label);
minTraingNum = 100000;
for i = minL : maxL
    Lindex = find(label==i);
    minNum = round(length(Lindex)*radio);
    if minNum < minTraingNum
        minTraingNum = minNum;
    end
end
display(sprintf('The min training num  is %d',minTraingNum));
trainIndex=[];
for i = minL : maxL
    Lindex = find(label==i);
    randomIndex = randperm(length(Lindex));
    Lindex = Lindex(randomIndex);
    if i == minL
        trainIndex = [trainIndex; Lindex(1:minTraingNum*4)];
    else
        trainIndex = [trainIndex; Lindex(1:minTraingNum)];
    end
end
testIndex = [1:1:length(label)];
testIndex(trainIndex)=[];
trainD = data(trainIndex);
trainL =  categorical(label(trainIndex));
testD =  data(testIndex);
testL =  categorical(label(testIndex));


% randomIndex = randperm( length(label));
% data = data(randomIndex);
% label = label(randomIndex);
% % 
% trainNum = round(radio *length(label));
% trainD = data(1:trainNum);
% trainL =  categorical(label(1:trainNum));
% testD =  data(trainNum:end);
% testL =  categorical(label(trainNum:end));



end

