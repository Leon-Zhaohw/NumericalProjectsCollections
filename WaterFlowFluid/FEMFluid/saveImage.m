function saveImage( frame, noFrame, outputFolder )
%SAVEIMAGE Save the given frame as an image in the output folder.
%   frame: matrix(yRes, xRes, colorChannel)
%   noFrame: simulation step the frame belongs to
%   outputFolder: the folder where the images are to be saved

    frame = squeeze(frame(:, :, :));

    if ~exist(outputFolder, 'dir')
        % Folder does not exist so create it.
        mkdir(outputFolder);
    end
    
    imwrite(frame, [outputFolder '\image' num2str(noFrame, '%04u') '.png']);

end

