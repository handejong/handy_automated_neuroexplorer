function cells = HAN_map_cells_to_atlas(source_file_path, varargin)
%MAP_CELLS_TO_ATLAS Maps inditified cells and maps the to the atlas.
%   Inputs are a picture (pixel map) of found cells produced by AIDA_histo
%   and the training_data MAT file which contains the tform object, which
%   detailes the transformation from the atlas slice to the microscopy
%   slice. We will not perform that transformation in reverse to map the
%   cells to the orginal ARA atlas coordinates.
%
%   This function will also add the original atlas slice as a Z coordinate
%   to each cell and label each cell with the brain region it was found.
%
%   Input arguments
%       - Path to the file created by AIDAhisto with all cells
%       - Patch to the .mat file with the training_data produced by
%         HAN_make reference
%       - Several optional input arguments:
%               - 'plot data': will plot the cells on top of the atlas
%               - 'save data': will save the data in a .txt file
%               - 'use cells': will instead of a source file use the cells
%               that are in the matrix one argument after this. Example:
%                   >> HAN_map_cells_to_atlas('1.jpg', 'use cells', cells);
%               were cells is a matrix [X, Y] with every row as a cell.
%
%
%   OUTPUT: an n by 6 matrix for n number of cells with the collowing
%   columns:
%       - X coordinate in atlas
%       - Y coordiante in atlas
%       - Z coordiante in atlas (Slice number)
%       - X coordiante in microscopy image
%       - Y coordiante in microscopy image
%       - Brain region identifyer (see ARA_legend for overview)
%
%   MAP_CELLS_TO_ATLAS is part of the Handy Automated Neuroexplorer (HAN)
%   made by Johannes de Jong at UC Berkeley. j.w.dejong@berkeley.edu

% Error handeling

% Deal with input arguments
plot_data = false;
save_data = false;
use_cells = false;
use_training_data = false;
force_un_edited = false;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case 'plot data' % plot the new data
            plot_data = true;
            
        case 'save data' % save the data
            save_data = true;
            
        case 'no editing'
            force_un_edited = true;
            
        case 'use cells'
            use_cells = true;
            old_cells = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'use training data'
            use_training_data = true;
            training_data = varargin{i+1};
            t_form = varargin{i+1}.tform;
            atlas = varargin{i+1}.atlas;
            varargin{i+1} = 'skipp';
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Get the paths
paths = HAN_get_paths(source_file_path);


% Open the file with the cells (either the edited one or the initial one)
% unless the user put in a matrix of cells themselves.
if ~use_cells
    
    if isfile(paths.cells_edited) && ~force_un_edited
        file_ID = fopen(paths.cells_edited);
    else
        file_ID = fopen(paths.cells);
    end


    % Read the cells from thi file
    cells_temp = textscan(file_ID, '%u %u', 'HeaderLines', 3);
    fclose(file_ID);
    old_cells(:,1) = cells_temp{1};
    old_cells(:,2) = cells_temp{2};
end


% Grab the t_form object, unless the user provided that
if ~use_training_data
    load(paths.training_data);
    t_form = training_data.tform;
    atlas = training_data.atlas;
end


% Grab the reference slice
if strcmp(atlas, '10um atlas')
    ref = imread(['Allen Reference Atlas/Slides 10um/slide_' num2str(training_data.ref_nr) '.tiff']);
elseif strcmp(atlas, '50um atlas')
    ref = imread(['Allen Reference Atlas/Slides 50um/slide_' num2str(training_data.ref_nr) '.tiff']);
end


% Note that the t_form object is just a transformation matrix from the
% reference to the source (microscopy image). The first column finds the
% new X coordiante as follows: new_X = old_X*t_form(1,1) +
% old_Y*t_form(2,1) + t_form(3,1). So the transformation is just a matrix
% multiplication: [old_x, old_y + 1] * t_form; (The 1 is for the
% translation).


% Add the ones used to calculate the bias (translation in this case)
old_cells(:,3) = ones(size(old_cells,1),1);


% Calculate the new cells
new_cells = double(old_cells)/t_form.T;


% If there is a scaling factor, apply that one as well
if isfield(training_data,'scale_factor')
    new_cells(:,1:2) = new_cells(:,1:2)./sqrt(training_data.scale_factor);
end


% Include the Z coordiante (slice number)
new_cells(:,3) = new_cells(:,3) * training_data.ref_nr;


% Include the old X and Y coordiantes
new_cells(:,4:5) = old_cells(:,1:2);


% Include an empty column where we will put the brain region
new_cells(:,6) = nan(size(new_cells, 1),1);


% For every cell find out in which brain region it is
for i=1:size(new_cells,1)
    X = round(new_cells(i,1));
    Y = round(new_cells(i,2));
    % Make sure the datapoint is actually in the atlas and not outside of
    % it (which is clearly wrong anyway)
    X = max([1, X]); X = min([size(ref,2), X]);
    Y = max([1, Y]); Y = min([size(ref,1), Y]);
    new_cells(i,6) = ref(Y, X);
end


% Save the data if requested
if save_data
    % get the file path
    filename = paths.cells_mapped;
    disp(['Saving data to: ' filename])    
    dlmwrite(filename, new_cells, 'delimiter', '\t');
end


% Visualize if requested
if plot_data
    figure; imagesc(ref); caxis([0 1000]); hold on; axis equal;
    scatter(new_cells(:,1), new_cells(:,2),'.','r');
    figure; image(training_data.source); axis equal;
end


% Put the data in the output arguemnt
cells = new_cells;


end

