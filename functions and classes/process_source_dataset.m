function flag = process_source_dataset(input_dataset, varargin)
%PROCESS_SOURCE_DATASET Takes care of preprocessing the input data
%   Input data to this function is a folder containing numbered image
%   files. The output data is a new folder containing the processed image
%   files. Process_source_dataset will open a GUI and present each image to
%   the user. The user accepts the image by pressing space or performs the
%   following opperations:
%
%   - f: will flip the image left-right
%
%   process_source_dataset is part of the Handy Automated Neuroexplorer
%   (HAN) made by Johannes de Jong at UC Berkeley. j.w.dejong@berkeley.edu

% Output flag, no succes
flag = 1;

% Error handeling 1, is it a folder
if ~isdir(input_dataset)
    error('Input should be a folder (include the path)')
end

% Find all the files in the folder
files = dir(input_dataset);

% Remove weird files that folders have no UNIX systems
mark_for_deletion = [];
for i=1:length(files)
    if strcmp(files(i).name, '..') || strcmp(files(i).name, '.') || strcmp(files(i).name, '.DS_Store')
        mark_for_deletion = [mark_for_deletion, i];
    end
end
files(mark_for_deletion) = [];

% Check if there are any other files in this folder then the source image
% files
% ...TODO...

% Figure out the path to the new folder and make it
new_dir = [input_dataset '_processed'];
mkdir(new_dir);

% Figure out the extention of the image files
index   = strfind(files(1).name,'.');
extension = files(1).name(index(end)+1:end);

% Main figure
main_figure = figure('KeyPressFcn',@key_press);
projection = image(imread([input_dataset '/1.' extension]));

% Variables
output_counter = 1;

% Loop through the images
for i = 1:length(files)
    c_filename = [num2str(i) '.' extension];
    c_filename_out = [num2str(output_counter) '.' extension];
    next_picture = false;
    
    % Try to import the next image
    try
        input_image = imread([input_dataset '/' c_filename]);
        projection.CData = input_image;
    catch
        warning([c_filename ' does not exist.'])
        next_picture = true;
    end
    
    % Wait untill we go to the next image
    while ~next_picture
        pause(0.1)
        if ~ishandle(main_figure)
            error('Figure deleted by the user')
        end
    end
    
end 

% OK we made it
flag = 1;
disp('Done')
delete(main_figure);


%%%%%%%%%%%%%%%%%%%%%% sub functions and callbacks %%%%%%%%%%%%%%%%%%%%%%%%

   function key_press(scr, ev)
        % Controls responses to key presses
        
        switch ev.Key
            
        
            case 'f' %Flip the image
                input_image = fliplr(input_image);
                projection.CData = input_image;
                
            case 'd' % Delete the image
                next_picture = true;
                
            case 'return' %Save the image
                imwrite(input_image, [new_dir '/' c_filename_out])
                next_picture = true;
                output_counter = output_counter + 1;
                
            % Unboud key
            otherwise
                disp(['Unbound key: ' ev.Key])
        end
        
        
    end





end
