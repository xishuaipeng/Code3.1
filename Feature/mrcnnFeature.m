function mrcnnFeature(videoPath,outPath)

bPath = pwd;
cd ('./Feature/')
command = sprintf('activate & python TriWapper.py -video_inpath "%s"  -txt_outpath "%s " & deactivate',videoPath,outPath);
system(command)
%py.importlib.import_module('D:\APP\Anaconda')
%pyversion D:\APP\Anaconda\\python.exe
% addpath ./
%py.importlib.import_module('skimage')
%py.TriWapper.loadMode()
%py.tss.pp(videoPath,outPath)
% commend = 'python TriWapper.py'
% system(commend)
cd (bPath)
end


