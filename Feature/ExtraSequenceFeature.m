function seqFeature = ExtraSequenceFeature(data, featureField, videoPath,startIndex, lastIndex,out_path )
    tableField = data.Properties.VariableNames;
    %%%%
    maxIndex = size(data,1);
    maxSeq = size(startIndex,1);
    indexUnique=[];
    for i=1:maxSeq
        indexUnique = [indexUnique,startIndex(i):lastIndex(i)];
    end
    indexUnique = unique(indexUnique);
    seqFeature = table(startIndex,lastIndex,'VariableNames',{'startIndex';'lastIndex'});
    %%%%%%%
    if sum(cellfun(@(x) strcmp({x},'VGGScene'),featureField)) > 0 & sum(cellfun(@(x) strcmp({x},'VGGScene'),tableField))== 0  
        frameUnique = data.frame(indexUnique);
        vgg_path = fullfile(out_path,'vggScene.mat');
        vggFeature = generate_vgg_scene(vgg_path,videoPath, indexUnique, frameUnique);
        seqFeature = feature2seq(vggFeature,seqFeature,'VGGScene');
    end
    if sum(cellfun(@(x) strcmp({x},'VGGObject'),featureField)) > 0 & sum(cellfun(@(x) strcmp({x},'VGGObject'),tableField))== 0  
        frameUnique = data.frame(indexUnique);
        vgg_path = fullfile(out_path,'vggObject.mat');
        vggFeature = generate_vgg_object(vgg_path,videoPath, indexUnique, frameUnique);
        seqFeature = feature2seq(vggFeature,seqFeature,'VGGObject');
    end
      
    if sum(cellfun(@(x) strcmp({x},'Heading'),featureField)) > 0 & sum(cellfun(@(x) strcmp({x},'Heading'),tableField))== 0  
%         save_path = fullfile(out_path,'Heading.mat');
        heading = data.GPS_heading(indexUnique);
        heading = [indexUnique', heading];
        seqFeature = feature2seq(heading,seqFeature,'Heading');
        seqFeature.Heading = cellfun(@(x) headingPro(x), seqFeature.Heading,'UniformOutput',false);
    end
    
    if sum(cellfun(@(x) strcmp({x},'Speed'),featureField)) > 0 & sum(cellfun(@(x) strcmp({x},'Speed'),tableField))== 0  
        Speed = data.speed(indexUnique);
        Speed = [indexUnique', Speed];
        seqFeature = feature2seq(Speed,seqFeature,'Speed');
        seqFeature.Speed = cellfun(@(x) x/100.0, seqFeature.Speed,'UniformOutput',false);
    end
    
    if sum(cellfun(@(x) strcmp({x},'Curvature'),featureField)) > 0 & sum(cellfun(@(x) strcmp({x},'Curvature'),tableField))== 0  
         
        x = data.GPS_lat(indexUnique);
        y = data.GPS_long(indexUnique);
        [x,y] = ll2utm(x,y);
        Curvature = [indexUnique',x,y];
        seqFeature = feature2seq(Curvature,seqFeature,'Curvature');
        seqFeature.Curvature = cellfun(@(x) curvaturePro(x), seqFeature.Curvature,'UniformOutput',false);
    end
    
  
%     if sum(cellfun(@(x) strcmp({x},'VGGObject'),featureField))>0 & sum(cellfun(@(x) strcmp({x},'VGGObject'),tableField))== 0  
%         vgg_path = fullfile(out_path,'vggObject.mat');
%         if ~exist(vgg_path)
%             curPath = mfilename('fullpath');
%             file  = strsplit( curPath,["/","\"]);
%             feature =0;
%             prototxtPath = fullfile(file{1:end-1}, 'VGG16', 'deploy_vgg16_places365.prototxt');
%             caffemodelPath = fullfile(file{1:end-1}, 'VGG16', 'vgg16_places365.caffemodel');
%             net = importCaffeNetwork(prototxtPath,caffemodelPath);
%             videoObj = VideoReader(videoPath);
%             vggFeature = zeros(langth_unique,4096+1);
%             vggFeature(:,1) = index_unique;
%             for i = 1 : langth_unique
%                 curFrame = index_unique(i);
%                 curImage = read(videoObj,curFrame);
%                 vggFeature(i,2:end) = vgg19Feature(curImage,net);
%             end
%             save(vgg_path, 'vggFeature');
%         else
%             clear vggFeature
%             load(vgg_path)
%         end
%         seq_index = vggFeature(:,1);
%         for i = 1:maxSeq
%             index = seq_index(i);
%             dura = lastIndex(i) - startIndex(i);
%             index = find(seq_index == index);
%             if(~seq_index( index + dura)== lastIndex(i))
%                 print('index is not matched in VGG mat!');
%                 return;
%             end
%             feature.vggScene{i}= vggFeature(index : index + dura, 2:end);
%         end

%     end
%     
%     
%     
%     
%     
%     if sum(cellfun(@(x) strcmp({x},'Curvature'),featureField))>0 
%        if sum(cellfun(@(x) strcmp({x},'Curvature'),tableField)) == 0
% 			[x,y] = ll2utm(data.GPS_lat,data.GPS_long);
%             kappa = Curvature(x,y);
%             kappa(kappa >2000) = 2000;
%             kappa(isnan(kappa)) = 2000;
%             kappa = kappa/2000;
%             kappa = smooth(kappa,0.2,'rloess');
%             %feature(:,index) = kappa;
%             data.Curvature =  kappa;
%        end
%     end
%     if sum(cellfun(@(x) strcmp({x},'Heading'),featureField))>0 
%        if sum(cellfun(@(x) strcmp({x},'Heading'),tableField))==0
%            smoothHeading = smooth(data.GPS_heading, 0.2, 'rloess');
%            cosHeading = cos(pi * smoothHeading/180);
%            cosHeading = (cosHeading+1)/2;
% %            gradHeading = smoothHeading(2:end) - smoothHeading(1:end-1);
% %            if any(gradHeading > 180)
% %                 smoothHeading = angRangeChange(smoothHeading);
% %            end
% %             meanH =   mean(smoothHeading);
% %            tHeading = smoothHeading - meanH ;
%            data.Heading =  cosHeading ;
%        end
%     end
%     if sum(cellfun(@(x) strcmp({x},'Speed'),featureField))>0
%         if sum(cellfun(@(x) strcmp({x},'Speed'),tableField))==0
%             speeding = data.speed/160;
%             data.Speed = speeding;
%         end
%          %feature(:,index) = speeding;
%     end
%     if sum(cellfun(@(x) strcmp({x},'Mask'),featureField))>0 
%        if  sum(cellfun(@(x) strcmp({x},'Mask'),tableField))==0
%             outPath = replace(videoPath,'.avi','_mask.txt');
%             outPath = [pwd,outPath(2:end)];
%             if ~isfile(outPath)
%                 tempVideo = [pwd,videoPath(2:end)];
%                 mrcnnFeature(tempVideo,outPath);
%             end       
%             fid = fopen(outPath);
%             maskFeature = [];
%             while ~feof(fid)
%                 tline = fgetl(fid);
%                 if isempty( tline)
%                     continue;
%                 end
%                 tline = split(tline,',');
%                 maskFeature = [maskFeature, tline];
%             end
%             fclose(fid);
%             maskFeature = cellfun(@(x) str2num(x),maskFeature);
%             maskFeature = maskFeature(2:end,:)';
%             %maskFeature = maskFeature(:,frameIndex)';
%             data.Mask = maskFeature;
%        end
%     end
    
end



function vggFeature = generate_vgg_object(feature_path, videoPath,indexUnique, frameIndex)
    if ~exist(feature_path)
                net = vgg19;
                n_case = length(frameIndex);
                videoObj = VideoReader(videoPath);
                vggFeature = zeros(n_case,4096+1);
                vggFeature(:,1) = indexUnique;
                for i = 1 : n_case
                    curFrame = frameIndex(i);
                    curImage = read(videoObj,curFrame);
                    vggFeature(i,2:end) = vgg19Feature(curImage,net);
                end
                save(feature_path, 'vggFeature');
     else
                clear vggFeature
                load(feature_path,'vggFeature')
    end
end

function vggFeature = generate_vgg_scene(feature_path, videoPath, indexUnique, frameIndex)
    if ~exist(feature_path)
                curPath = mfilename('fullpath');
                file  = strsplit( curPath,["/","\"]);
                prototxtPath = fullfile(file{1:end-1}, 'VGG16', 'deploy_vgg16_places365.prototxt');
                caffemodelPath = fullfile(file{1:end-1}, 'VGG16', 'vgg16_places365.caffemodel');
                net = importCaffeNetwork(prototxtPath,caffemodelPath);
                n_case = length(frameIndex);
                videoObj = VideoReader(videoPath);
                vggFeature = zeros(n_case,4096+1);
                vggFeature(:,1) = indexUnique;
                for i = 1 : n_case
                    curFrame = frameIndex(i);
                    curImage = read(videoObj,curFrame);
                    vggFeature(i,2:end) = vgg16Feature(curImage,net);
                end
                save(feature_path, 'vggFeature');
     else
                clear vggFeature
                load(feature_path,'vggFeature');
    end
end

function seqFeature = feature2seq(modelFeature,seqFeature,field)
        modelIndex = modelFeature(:,1);
        maxSeq = size(seqFeature,1);
        for i = 1:maxSeq 
            seq_start = seqFeature.startIndex(i);
            seq_duration = seqFeature.lastIndex(i)- seq_start;
            index = find( modelIndex == seqFeature.startIndex(i));
            if(length( modelIndex( index + seq_duration) == seqFeature.lastIndex(i))==0)
                fprintf('index is not matched in %s mat!',field );
                return;
            end
        seqFeature.(field){i}= modelFeature(index : index + seq_duration, 2:end);
        end
end

function headFeature = headingPro(x)
    cosx = cos(x*pi/180);%[-1,1]
    meanHeading = mean(cosx);%[-1,1]
    cosx = cosx - meanHeading;%[-2,2]
    cosx = abs(cosx/2);%[0,1]
    headFeature = cosx;
end


function curvatureFeature = curvaturePro(data)
    x = data(:,1);
    y = data(:,2);
    x_ = diff(x,1);
    x_ = x_(2:end);
    x__ = diff(x,2);
    X = x(3:end);
    y_ = diff(y,1);
    y_ = y_(2:end);
    y__ = diff(y,2);
    Y = y(3:end);  
    curvatureFeature = abs(x_.* y__ - x__ .* y_)./(x_.^2 + y_.^2).^(3/2);
    curvatureFeature = [curvatureFeature(1);curvatureFeature(1);curvatureFeature];
    curvatureFeature(isnan(curvatureFeature))=2000;
    curvatureFeature(curvatureFeature>2000)=2000;
    curvatureFeature=curvatureFeature/2000;

end






function angleList = angRangeChange(angleList, varargin)
% change angle range
% angleList = angRangeChange(angleList, outMax)
% angleList = angRangeChange(angleList)
% written by Ruirui Liu
if nargin == 2
    outMax = varargin{1};
for i = 1: length(angleList)
    if angleList(i)<outMax-360
        n = ceil((outMax-360-angleList(i))/360);
        angleList(i) = angleList(i)+n*360;
    elseif angleList(i)>outMax
        n = ceil((angleList(i)-outMax)/360);
        angleList(i) = angleList(i)-n*360;
    end
end
elseif nargin == 1
    for i = 2:length(angleList)
        if abs(angleList(i) - angleList(i-1)) > 180
            angleList(i) = angleList(i) - ...
                360*round((angleList(i) - angleList(i-1))/360);
        end
%         theta1 = angleList(i-1)/180*pi;
%         theta2 = angleList(i)/180*pi;
%         angleList(i) = angleList(i-1) + ...
%             asin(sin(theta2-theta1))/pi*180;
    end
else 
end

end % EOF

