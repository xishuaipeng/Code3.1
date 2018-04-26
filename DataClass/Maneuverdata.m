classdef Maneuverdata < Dataset       

    methods(Static)
        function y = LabelManeuver(eventFrame, x)
            [r,c] = size(x);
            frame = zeros(r,2);
            for i = 1:r
               frame(i,:) = [min(x{i,1}(:,1)), max(x{i,1}(:,1))]; 
            end
            y = zeros(r,1); 
            eventNum = size(eventFrame,1);
            for j = 1:eventNum
                  event_start = eventFrame(j,1);
                  event_last =  eventFrame(j,2);
                  event = find(eventFrame(j,3:end)==1);
                for i= 1:r
                    seq_start = min(x{i}(:,1));
                    seq_last = max(x{i}(:,1));     
                    start = max(event_start, seq_start);
                    last = min(event_last, seq_last);
                    duration = last - start;
                    radio_event = duration /(event_last - event_start);
                    radio_seq  =  duration /(seq_last - seq_start);
                    radio = max(radio_event,radio_seq);
                    if radio > 0.7
                        y(i) = event;
                    end
                end
            end   
        end

        function encodeevent(fid, cur_dir, score,descript)
             fwrite(fid, sprintf('%s,%f,%s \n',cur_dir,score,descript));  
        end
        
        function [cur_dir, score,descript] = decodeevent(line_ex)
            component = split(line_ex,',');
            cur_dir = component{1};
            score = str2num(component{2});
            descript = component{3};
        end

    end
    methods
        function obj = Maneuverdata(varargin)   
            obj@Dataset()
            obj.localParameters.addParamValue('seq_window', 0.05);
            obj.localParameters.addParamValue('seq_step', 0.01);
            obj.localParameters.addParamValue('gene_feature', {'Curvature','VGGObject','Heading','Speed'});
            obj.localParameters.addParamValue('load_feature', {'Curvature','VGGObject','Heading','Speed'});
            obj.localParameters.addParamValue('start_padding', 0.1);
            obj.localParameters.addParamValue('last_padding', 0.02);
            obj.localParameters.addParamValue('out_dir', '');
            obj.localParameters.addParamValue('syn_file_path','');
            obj.localParameters.addParamValue('sample_data_path', '');
            obj.localParameters.addParamValue('sample_event_path','');
            obj.localParameters.addParamValue('event_checklist_path','');
            obj.localParameters.addParamValue('valid_sample_event_path','');  
            obj.localParameters.addParamValue('event_data_dir',''); 
            obj.localParameters.addParamValue('seq_y_path',''); 
            obj.localParameters.addParamValue('seq_x_path',''); 

            obj.Initialization(varargin{:});          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%update parameters
            input_path = obj.localParameters.Results.input_path;
            data_id = obj.localParameters.Results.data_id;
            if isempty(obj.localParameters.Results.out_dir);obj.localParameters = obj.update...
                    (obj.localParameters,'out_dir', sprintf('%s/processed_data',input_path));end
            mkdir_if_not_exist(fullfile(obj.localParameters.Results.out_dir,data_id));
            if isempty(obj.localParameters.Results.syn_file_path);obj.localParameters = obj.update...
                    (obj.localParameters,'syn_file_path', sprintf('%s/sycnFile.txt',input_path));end

            out_dir = obj.localParameters.Results.out_dir;
            samplingProperty = obj.localParameters.Results.samplingProperty;
            samplingStep = obj.localParameters.Results.samplingStep;
            if isempty(obj.localParameters.Results.sample_data_path);obj.localParameters = obj.update...
                    (obj.localParameters,'sample_data_path', sprintf('%s/%s/sampleData_%s_%.3f.mat',out_dir,data_id,samplingProperty,samplingStep));end
            if isempty(obj.localParameters.Results.sample_event_path);obj.localParameters = obj.update...
                    (obj.localParameters,'sample_event_path', sprintf('%s/%s/sampleEvent_%s_%.3f.mat',out_dir,data_id,samplingProperty,samplingStep));end
            if isempty(obj.localParameters.Results.valid_sample_event_path);obj.localParameters = obj.update...
                    (obj.localParameters,'valid_sample_event_path', sprintf('%s/%s/ValidsampleEvent_%s_%.3f.mat',out_dir,data_id,samplingProperty,samplingStep));end
            if isempty(obj.localParameters.Results.event_checklist_path);obj.localParameters = obj.update...
                    (obj.localParameters,'event_checklist_path', sprintf('%s/%s/%s.list',out_dir,data_id,data_id));end
            
            if isempty(obj.localParameters.Results.event_data_dir);obj.localParameters = obj.update...
                (obj.localParameters,'event_data_dir', sprintf('%s/%s/Event',out_dir,data_id));end
            seq_window = obj.localParameters.Results.seq_window;
            seq_step = obj.localParameters.Results.seq_step;
            if isempty(obj.localParameters.Results.seq_y_path);obj.localParameters = obj.update...
                (obj.localParameters,'seq_y_path', sprintf('%s/%s/seqY_%f_%f.mat',out_dir,data_id,seq_window,seq_step));end     
            if isempty(obj.localParameters.Results.seq_x_path);obj.localParameters = obj.update...
                (obj.localParameters,'seq_x_path', sprintf('%s/%s/seqX_%f_%f.mat',out_dir,data_id,seq_window,seq_step));end     
       
        end
        function [obj, sample_data] = resampLogdata(obj)
            sample_data_path = obj.localParameters.Results.sample_data_path;
            videoPath = obj.localParameters.Results.videoPath;
            vidObj = VideoReader(videoPath);
            obj.localParameters = obj.update(obj.localParameters,'frameRate',vidObj.FrameRate); 
            maxFrame = vidObj.NumberOfFrames;
            [obj,logData] = obj.readLogdata();
            syn_file_path = obj.localParameters.Results.syn_file_path;
            data_id = obj.localParameters.Results.data_id;
            timeDelayforVideo = obj.resync(data_id, logData, syn_file_path);
            obj.localParameters = obj.update(obj.localParameters,'timeDelayforVideo',timeDelayforVideo);     
            if exist(sample_data_path,'file')
                clear sample_data
                load(sample_data_path,'sample_data');
            else
                [obj,sample_data] = resampLogdata@Dataset(obj, logData);
                sample_data(sample_data.frame>maxFrame,:)=[];
                sample_data(sample_data.frame<1,:)=[];
                
                save(sample_data_path,'sample_data' );
            end
        end      
        function sample_event = resample_event_from_data(obj, data)
            sample_event_path = obj.localParameters.Results.sample_event_path;
            event_field = obj.localParameters.Results.event_field;
            if exist(sample_event_path,'file')
                clear sample_event
                load(sample_event_path,'sample_event');
            else
                %function part 
                sample_event= obj.readEvent();
                num_event = size(sample_event,1);
                dataframe = data.frame;
                maxFrame =  max(dataframe);
                for eventNum = 1: num_event  
                    startFrame = sample_event.StartTime(eventNum);
                    endFrame = sample_event.EndTime(eventNum);
                    eventIndex = find(sample_event{eventNum,event_field} == 1);
%                     eventIndex = find(event_frame(eventNum,event_field) == 1);
                    compIndex = min(1,dataframe - startFrame);
                    compIndex(compIndex== 1) = - maxFrame;
                    [~, startIdx_] = max(compIndex);
                    %end index
                    compIndex = max(-1,dataframe - endFrame);
                    compIndex(compIndex == -1) =  maxFrame;
                    [~, endIdx_] = min(compIndex);
                    sample_event.StartsampleIndex(eventNum) = startIdx_;
                    sample_event.LastsampleIndex(eventNum) = endIdx_; 
                    sample_event.StartsampleFrme(eventNum) = dataframe(startIdx_);
                    sample_event.LastsampleFrame(eventNum) = dataframe(endIdx_); 
                    sample_event.Durationtime(eventNum)  = data.time(endIdx_) - data.time(startIdx_);
                    sample_event.Durationdistance(eventNum)  =  data.distance(endIdx_) - data.distance(startIdx_);
                    sample_event.Type(eventNum) = eventIndex;
                    sample_event.Index(eventNum) = eventNum;
                end  
                save(sample_event_path,'sample_event');
            end
            
            event_checklist_path = obj.localParameters.Results.event_checklist_path;
            if ~exist(event_checklist_path,'file')
                list_writer = fopen(event_checklist_path,'w+');
                try
                num_event = size(sample_event,1);
                for eventNum = 1: num_event 
                    eventIndex = sample_event.Type(eventNum);
                    cur_dir = sprintf('%d_%s',eventNum, event_field{eventIndex});
                    obj.encodeevent(list_writer, cur_dir, 1,'NEEDCHECK');  
                end
                catch
                    fclose(list_writer);
                    display('write list file wrong!')
                    return; 
                end
                fclose(list_writer);
            end
        end  
        function sample_event = load_resample_event(obj,sample_data)
            valid_sample_event_path = obj.localParameters.Results.valid_sample_event_path;
            if exist(valid_sample_event_path)
                clear sample_event;
                load(valid_sample_event_path,'sample_event');
            else
                event_checklist_path = obj.localParameters.Results.event_checklist_path;
                sample_event_path = obj.localParameters.Results.sample_event_path;
                if ~exist(event_checklist_path,'file')|~exist(sample_event_path,'file')
                    sample_event = obj.resample_event_from_data(sample_data);
                else
                    clear sample_event
                    load(sample_event_path,'sample_event'); 
                end
                valid_reader = fopen(event_checklist_path,'r');
                index = 0;
                 while ~feof(valid_reader)
                        tline = fgetl(valid_reader);
                        if isempty(tline);continue;end
                        index = index +1;
                        [dirStr,score,description] = obj.decodeevent(tline);
                        if score < 0.5
                            sample_event(index,:) = [];
                            continue;
                        else
                            dirStr = replace(dirStr,'Event\','');
                            if strcmp(description,'NEEDCHECK');display([dirStr,' is not chenked']);end  
                        end
                 end
                 save(valid_sample_event_path,'sample_event');   
                 fclose(valid_reader);
            end
        end
        function sequenceZ = sequence_from_event(obj, sample_data)
            seq_y_path =  obj.localParameters.Results.seq_y_path;
            if exist(seq_y_path,'file')
                clear sequenceZ
                load(seq_y_path,'sequenceZ');
             else
                event_frame = obj.load_resample_event(sample_data);
                num_event = size(event_frame,1);
                start_padding = obj.localParameters.Results.start_padding;
                last_padding = obj.localParameters.Results.last_padding;
                samplingStep = obj.localParameters.Results.samplingStep;
                start_padding = start_padding/samplingStep;
                last_padding = last_padding/samplingStep;

                seq_window = obj.localParameters.Results.seq_window;
                seq_step = obj.localParameters.Results.seq_step;
                samplingStep = obj.localParameters.Results.samplingStep;
                seq_window = floor(seq_window/samplingStep);
                seq_step = floor(seq_step/samplingStep);
                maxIndex = size(sample_data,1);
                sequence_index = [];
                for i = 1 : num_event
                    startIndex = event_frame.StartsampleIndex(i);
                    lastIndex =  event_frame.LastsampleIndex(i);
                    eventIndex = event_frame.Index(i);
                    eventLabel = event_frame.Type(i);
                    startPadding = max(1,round(startIndex - start_padding));
                    lastPadding = min(maxIndex,round(lastIndex + last_padding));
                    end_start_index = (lastPadding - seq_window);
                    sequence_start = [startPadding:seq_step:end_start_index];
                    if isempty(sequence_start)
                        fprintf('Event %d is too short,please recheck!',eventIndex);
                        continue;
                    end
                    if ~(sequence_start(end)==end_start_index)
                        sequence_start = [sequence_start,end_start_index];
                    end
                    for j = 1: length(sequence_start)
                        seq_start = sequence_start(j);
                        seq_end = sequence_start(j) + seq_window ;
%                         bestIOU = 0 ;
%                         bestIndex = 0;
%                         for k = 1: num_event
%                             overlap = max(0, min(seq_end, event_frame.LastsampleIndex(k))  - max(seq_start, event_frame.StartsampleIndex(k)));
%                             if overlap>bestIOU;bestIOU = overlap;end  
%                         end
%                         if bestIndex
                        SE2EE_time =  sample_data.time(lastIndex) - sample_data.time(seq_end);
                        SE2EE_dis = (sample_data.distance(lastIndex) - sample_data.distance(seq_end))*1000 ;
                        SE2EB_time =  sample_data.time(startIndex) - sample_data.time(seq_end);
                        SE2EB_dis = (sample_data.distance(startIndex) - sample_data.distance(seq_end))*1000 ;
                        overlap = max(0, min(seq_end, lastIndex)  - max(seq_start, startIndex));
                        SinE =  overlap/(lastIndex -startIndex );
                        EinS =  overlap/(seq_end -seq_start );
                        seqLabel = eventLabel;
                        if SinE <0.3 & EinS<0.3 ;seqLabel = 0; end
                        sequence_index = [sequence_index; eventIndex,seqLabel,seq_start,seq_end, SE2EE_time,SE2EE_dis,SE2EB_time,SE2EB_dis ,SinE,EinS];
                    end
                end
                sequenceZ = array2table(sequence_index);
                sequenceZ.Properties.VariableNames = {'eventIndex','eventType','startIndex','lastIndex','seqEnd2eventEnd_time','seqEnd2eventEnd_distance', 'seqEnd2eventBeg_time','seqEnd2eventBeg_distance','SinE','EinS'};  
                save(seq_y_path, 'sequenceZ');
            end
         end
        function WriteEvent(obj)
            [obj, sample_data] = obj.resampLogdata();
            start_padding = obj.localParameters.Results.start_padding;
            last_padding = obj.localParameters.Results.last_padding;
            samplingStep = obj.localParameters.Results.samplingStep;
            start_padding = start_padding/samplingStep;
            last_padding = last_padding/samplingStep;
            event_checklist_path = obj.localParameters.Results.event_checklist_path;
            valid_reader = fopen(event_checklist_path,'r');
            videoPath = obj.localParameters.Results.videoPath;
            event_data_dir = obj.localParameters.Results.event_data_dir;
            vidObj = VideoReader(videoPath);
            frameRate = obj.localParameters.Results.frameRate;
            sample_event = obj.load_resample_event(sample_data);
            index = 0;
             while ~feof(valid_reader)
                    tline = fgetl(valid_reader);
                    if isempty(tline);continue;end
                    index = index +1;
                    [dirStr,score,description] = obj.decodeevent(tline);
                    if score < 0.5;continue; end
                    dirStr = replace(dirStr,'Event\','');
                    if strcmp(description,'NEEDCHECK');display([dirStr,' is not chenked']);end
                    startFrame = sample_event.StartsampleIndex(index);
                    endFrame = sample_event.LastsampleIndex(index);
                    startIdx = max(1,round(startFrame - start_padding));
                    endIdx = min(length(sample_data.frame), round(endFrame + last_padding));
                    if startIdx >= endIdx
                        continue;
                    end
                    segData = sample_data(startIdx : endIdx, :);
                    savePath = fullfile( event_data_dir,dirStr);
                    mkdir_if_not_exist(savePath);
                    matPath = fullfile(savePath, sprintf('%d_%d_%.1f.mat',start_padding,last_padding, 1000*samplingStep));
                    save(matPath, 'segData');  
                    %write
                    segVidObj = VideoWriter(fullfile(savePath , sprintf('%d_%d_%.1f.avi',start_padding,last_padding, 1000*samplingStep)));
                    segVidObj.FrameRate = frameRate;
                    open(segVidObj);
                    segFrame= segData.frame;
                    % for experiment
                    for frameIndex = 1 : length(segFrame) 
                        writeVideo(segVidObj,read(vidObj,segFrame(frameIndex)));
                    end
                    close(segVidObj);
             end
             fclose(valid_reader);
        end
        function [seqX,sequenceZ] = EventSequence(obj)
            seq_x_path = obj.localParameters.Results.seq_x_path;

            [obj, sample_data] = obj.resampLogdata();
            data_id =  obj.localParameters.Results.data_id;
            out_dir = obj.localParameters.Results.out_dir;
            gene_feature = obj.localParameters.Results.gene_feature;
            load_feature = obj.localParameters.Results.load_feature;
            videoPath = obj.localParameters.Results.videoPath;
            sequenceZ = obj.sequence_from_event(sample_data);
%             if sequence
            if exist(seq_x_path)
                clear seqFeature;
                load(seq_x_path,'seqFeature');
            else
                seqFeature = ExtraSequenceFeature(sample_data, gene_feature, videoPath,sequenceZ.startIndex, sequenceZ.lastIndex,fullfile(out_dir,data_id));
                save(seq_x_path,'seqFeature');
            end
            seqFeature = table2array( seqFeature(:,load_feature));
            nCase = size(seqFeature,1);
            seqX = cell(nCase,1);
            for i = 1:nCase; seqX(i,1)={cat(2,seqFeature{i,:})}; end
            seqX = cellfun(@(x) x', seqX,'UniformOutput',false);
        end
    end
    
end