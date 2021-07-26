function saveVideo( visualizationSettings, deleteImages )
%SAVEVIDEO Create a video out of all .png images in the output folder.
%   visualisationSettings: struct with these fields:
%       outputFolder: folder where the resulting video should be saved
%       fps: frames per second of the video
%       outputFile: the name of the resulting video
%   deleteImages: bool delete or keep source images

    display('Saving video...');
    
    [~,~,~] = mkdir(visualizationSettings.outputFolder);
    
    if (exist([visualizationSettings.outputFolder '\' visualizationSettings.outputFile], 'file'))
        delete([visualizationSettings.outputFolder '\' visualizationSettings.outputFile]);
    end
    
    videoFile = [visualizationSettings.outputFolder '\' visualizationSettings.outputFile];
    
    % initialize video writer
    writerObj = VideoWriter(videoFile);
    writerObj.FrameRate = visualizationSettings.fps;
    open(writerObj);
    
    % open every single output frame and add it to the video
    file_list_frames = dir([visualizationSettings.outputFolder '/image*.png']);
    for j = 1:numel(file_list_frames)
        frame = imread([visualizationSettings.outputFolder '/' file_list_frames(j).name]);
        writeVideo(writerObj, frame);
    end
    
    if deleteImages
        delete([visualizationSettings.outputFolder '\image*.png']);
    end
    
    display('Finished!');
    
end

