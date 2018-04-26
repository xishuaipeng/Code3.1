classdef  Dataset
    
properties
localParameters;
logTable;
logData;
sampleData;
eventLabel;
eventFrame;
end

methods(Static)
  
    function parameter = update(parameter,varargin)
      n_case = length(varargin);
     
      if n_case == 0
          return ;
      end
      field = parameter.Parameters;
      key_index = [1:2:n_case];
      value_index = [2:2:n_case];
      mlist = parameter.Results;
      n_pair = length(key_index);
      for i= 1 : n_pair
          try 
              key = varargin{key_index(i)};
              index = cellfun(@(x) strcmp(x, key), field);
              index = sum(index);
              if index>0
                mlist.(key) = varargin{value_index(i)};    
              end
          catch
              continue;   
          end   
      end
      parse(parameter,mlist);
    end  
    function sequence = table2sequence(table, window, step)
        [r,c] = size(table);
        numSequence = floor((r-window)/step)+1;
        sequence = cell(numSequence,1);    
       for i= 1:numSequence
            start = (i-1)*step;
            last = start + window;
            seq = table(start+1:last, :);
            sequence{i} = seq;
       end
    end
    function timeDelayforVideo = resync(data_id, lod_data, file_path)
            timeDelayforVideo = checkIfSync(lod_data, data_id, file_path);
      end
    function  eventFrame= readeventlabel(eventField, labelPath,frameRate)
          eventLabel  = readtable(labelPath);
          num_case = size(eventLabel,1);
%           field = {'StartTime','EndTime', eventField{:}};
          eventFrame = zeros(num_case, 2+length(eventField));
          for i=1: num_case
            try
                startIndex = round(timeParser(eventLabel.StartTime{i})* frameRate);
                endIndex =  round(timeParser(eventLabel.EndTime{i})* frameRate);
            catch 
                 startIndex = round(seconds(eventLabel.StartTime(i))* frameRate);
                 endIndex =  round(seconds(eventLabel.EndTime(i))* frameRate); 
            end
            event= eventLabel{i,eventField};
            if sum(event)>0
                eventFrame(i,:) = [startIndex,endIndex, event ];
            end
          end
          eventFrame(sum(eventFrame,2)==0,:)=[];
          eventFrame = array2table(eventFrame);
          eventFrame.Properties.VariableNames ={'StartTime','EndTime', eventField{:}};
    end  
 end
methods
    function obj = Dataset()
        obj.localParameters=inputParser;
        obj.localParameters.KeepUnmatched = true;  
    end
    function obj = Initialization(obj, varargin)
        obj.localParameters.addParamValue('data_id', '');
        obj.localParameters.addParamValue('input_path',  '');
        obj.localParameters.addParamValue('log_field', {'time','speed','GPS_long','GPS_lat','GPS_heading','distance'});
        obj.localParameters.addParamValue('event_field', {'TurnLeft','TurnRight','LaneChangeLeft','LaneChangeRight'});
        obj.localParameters.addParamValue('timeDelayforVideo',0);
        obj.localParameters.addParamValue('samplingProperty','distance');
        obj.localParameters.addParamValue('samplingStep',0.002);
        obj.localParameters.addParamValue('logRate',0);
        obj.localParameters.addParamValue('logPath', '');
        obj.localParameters.addParamValue('labelPath','');
        obj.localParameters.addParamValue('videoPath','');
        obj.localParameters.addParamValue('frameRate',0);
        obj.localParameters.parse(varargin{:}); 
        %update
        dataid = obj.localParameters.Results.data_id;
        datapath = obj.localParameters.Results.input_path;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isempty(obj.localParameters.Results.logPath);obj.localParameters = obj.update(obj.localParameters,...
                'logPath', sprintf('%s/raw_data/%s/%s_datalog.Csv',datapath,dataid,dataid));end
        
        if isempty(obj.localParameters.Results.labelPath);obj.localParameters =obj.update(obj.localParameters,...
                'labelPath', sprintf('%s/raw_data/%s/%s_labelresult.Csv',datapath,dataid,dataid));end
        
        if isempty(obj.localParameters.Results.videoPath);obj.localParameters = obj.update(obj.localParameters,...
        'videoPath', sprintf('%s/raw_data/%s/%s_video.avi',datapath,dataid,dataid));end
    end 
    function [obj, logData] = readLogdata(obj)
        logPath = obj.localParameters.Results.logPath;
        logField = obj.localParameters.Results.log_field;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pathout = [logPath(1:end-4),'.mat'];
        logTable= TRI_GPS_extract_oneTrip(logPath, pathout);
%         logTable  = TRI_GPS_extract_oneTrip(logPath);
        logTable = logTable(:,logField);
        %table 2 data
        logData  = table();
        for i = 1:size(logTable ,2)
            logData (:,i) = table(str2num(char(table2array(logTable (:,i)))));
        end
        logData.Properties =logTable.Properties;
        
        logRate = logData.time(2) - logData.time(1) ;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        obj.localParameters = obj.update(obj.localParameters,'logRate', logRate);
    end  
    function [obj,sampleData] = resampLogdata(obj, logData)
        frameRate = obj.localParameters.Results.frameRate;
        timeDelayforVideo = obj.localParameters.Results.timeDelayforVideo;
        samplingStep = obj.localParameters.Results.samplingStep;
        samplingProperty = obj.localParameters.Results.samplingProperty;   
        if isempty(logData)
            print('logData is empty!')
        end
        [num_case, num_feature] = size(logData);
        segVector = logData.(samplingProperty);
        max_sample = 1 + floor((segVector(end) - segVector(1)) / samplingStep );
        %          %find the time field index
        timeIndex = cellfun(@(x) strcmp(x,'time'),logData.Properties.VariableNames);
        timeIndex = find(timeIndex==1);
        %         %construct the array
        tempData = zeros(max_sample, num_feature);
        frameIndex = zeros(max_sample, 1);
        interValue = [segVector(1) : samplingStep : segVector(end)];
        tempData(1,:) = logData{1,:};
        frameIndex(1)=  round((tempData(1,timeIndex)  - timeDelayforVideo )* frameRate )+1;
        index = 1;  
        for  i=2:num_case
            current =interValue(index) + samplingStep;%tempData(index, segIndex) 
            if logData.(samplingProperty )(i) >= current
                x0 = logData.( samplingProperty )(i-1);
                x1 = logData.( samplingProperty )(i);
                x = current;
                theta = (x-x0)/(x1-x0);
                index = index + 1; 
                tempData(index,1:end) = logData{i-1,1:end}+ theta * (logData{i,1:end}-logData{i-1,1:end});  
                %transfer time to frame index
                frameIndex(index)=  round((tempData(index,timeIndex)  - timeDelayforVideo )* frameRate )+1;
            end   
        end  
        sampleData = array2table(tempData);
        sampleData.Properties = logData.Properties;
        sampleData.frame = frameIndex;
    
    end

    function eventFrame = readEvent(obj)
        event_field = obj.localParameters.Results.event_field;
        frameRate = obj.localParameters.Results.frameRate;
        labelPath = obj.localParameters.Results.labelPath;
        eventLabel  = readtable(labelPath);
        num_case = size(eventLabel,1);
%           field = {'StartTime','EndTime', eventField{:}};
          eventFrame = zeros(num_case, 2+length(event_field));
          for i=1: num_case
            try
                startIndex = round(timeParser(eventLabel.StartTime{i})* frameRate);
                endIndex =  round(timeParser(eventLabel.EndTime{i})* frameRate);
            catch 
                 startIndex = round(seconds(eventLabel.StartTime(i))* frameRate);
                 endIndex =  round(seconds(eventLabel.EndTime(i))* frameRate); 
            end
            event= eventLabel{i,event_field};
            if sum(event)>0
                eventFrame(i,:) = [startIndex,endIndex, event ];
            end
          end
          eventFrame(sum(eventFrame,2)==0,:)=[];
          eventFrame = array2table(eventFrame);
          eventFrame.Properties.VariableNames ={'StartTime','EndTime', event_field{:}};
    end  
    function obj = labelSampledata(varargin)
          obj = varargin{1};
          if isempty(obj.sampleData)
              'no sampleData please generate it first'   
          end
          if  isempty(obj.eventLabel)
              obj = obj.readEventlabel();
              
          end
          maxFrame = obj.sampleData.frame(end);
          tempFLabel = zeros(maxFrame,1);
          num_case = size(obj.eventFrame,1);
%           eventFrame = zeros(num_case, 2+length(eventField));
          for i=1: num_case
            startIndex = obj.eventFrame(i,1);
            endIndex =  obj.eventFrame(i,2);
            event= obj.eventFrame(i,3:end);
            event = find(event>0);
            tempFLabel(startIndex:endIndex,:) =  repmat(event,endIndex - startIndex+1,1);
          end
          obj.sampleData.label = tempFLabel(obj.sampleData.frame,:);
    end 
    function obj = reSync(varargin)
        obj = varargin{1};
        if numel(varargin)==1
            obj.timeDelayforVideo = checkIfSync(obj.logData);
        elseif numel(varargin) == 4
            obj.timeDelayforVideo = checkIfSync(obj.logData,varargin{2},varargin{3},varargin{4});
        elseif numel(varargin) == 2
            obj.timeDelayforVideo = checkIfSync(obj.logData, obj.dataID, varargin{2});
        end
     end
    function obj = checkLabel(varargin)
      obj = varargin{1};
      x = obj.sampleData.GPS_long;
      y = obj.sampleData.GPS_lat;
      maxLength = length(x);
      label = obj.sampleData.label+1;
      %label =label - min(label);
      classes = max(label);
      radio = jet( classes);
      figure; hold on;
      for i = 1: classes
          index = find(label == i);
          p_x = x(index);
          p_y = y(index);
          color = radio(i,:);
          scatter(p_x,p_y, 'MarkerFaceColor',color,'MarkerEdgeColor',color);
          
      end
     legend(['No event',obj.eventField]);
      
      
      
%       maxLength = length(x);
%       event_start = round((obj.timeParser(obj.eventLabel.StartTime )+ obj.timeDelayforVideo ) * logRate)+1;
%       event_end =  round((obj.timeParser(obj.eventLabel.EndTime) + obj.timeDelayforVideo )* logRate)+1;
    %           event_time =[ event_start, event_end];
%       eventNum = length(event_start);
%       figure();hold on;
%       plot(x,y,'.k');
%       colorBar = ['r', 'g', 'y', 'b'];
%       event = table2array(obj.eventLabel(:,obj.eventField));
%       
%       
%       
%        for i =1 : eventNum
%            hold on;
%            label_index = event(i,:)
%            label_index = find(label_index==1);
%            if ~isempty(label_index)
%                obj.frameLabel(event_start(i):event_end(i))= obj.eventField{label_index};
%                plot(x(event_start(i):event_end(i)),y(event_start(i):event_end(i)),'.', ...
%         'Color', colorBar(label_index), 'MarkerSize', 6);
%            end
%        end               
    end  
end

end