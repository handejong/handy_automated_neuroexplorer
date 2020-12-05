function flag = HAN_update_atlas_data(source_image_path, training_data)
%HAN_UPDATE_ATLAS_DATA updates a .mat file with info about which atlas
%slide is linked to which brain slide.
%   ... TODO 

% Output flag
flag = 1;

% Grab the path
[path, filename, ext] = fileparts(source_image_path);
atlas_data_path = [path '/' 'atlas_data.mat'];

% Load if exist, otherwise make a new one
if isfile(atlas_data_path)
    load(atlas_data_path,'atlas_data');
    
    % Check if there is allready info about this slide
    found_slice = false;
    for i=1:length(atlas_data)
        if strcmp(atlas_data(i).source, [filename ext])
            
            % Error handeling (found slice twice)
            if found_slice
                warning('Two coppies of this slice in the atlas data file.')
                atlas_data(i) = [];
            else
                found_slice = true;
                atlas_data(i).slide_number = str2num(filename);
                atlas_data(i).source = [filename ext];
                atlas_data(i).atlas_type = training_data.atlas;
                atlas_data(i).ref_nr = training_data.ref_nr;
                atlas_data(i).tform = training_data.tform;
            end
            
        end
    end
        
    % If we did not find the slice
    if ~found_slice
        atlas_data(end+1).slide_number = str2num(filename);
        atlas_data(end).source = [filename ext];
        atlas_data(end).atlas_type = training_data.atlas;
        atlas_data(end).ref_nr = training_data.ref_nr;
        atlas_data(end).tform = training_data.tform;
    end
else
    atlas_data.slide_number = str2num(filename);
    atlas_data.source = [filename ext];
    atlas_data.atlas_type = training_data.atlas;
    atlas_data.ref_nr = training_data.ref_nr;
    atlas_data.tform = training_data.tform;
end

% Save the new atlas_data
save(atlas_data_path, 'atlas_data');
disp('Atlas data updated')

% Done
flag = 0;
