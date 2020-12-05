function flag = HAN_show_on_template(whole_brain_file, atlas_type, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Output flag
flag = 1;

% Deal with input arguments
plot_data = false;
save_data = false;
level = 1;
min_cells = 1;
start_slice = 100;
orientation = 'coronal';
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case 'level'
            level = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'min cells'
            min_cells = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'start slice'
            start_slice = varargin{i+1};
            vargin{i+1} = 'skipp';
            
        case 'sagittal'
            orientation = 'sagittal';
            
        case 'coronal'
            orientation = 'coronal';
            
        case 'horizontal'
            orientation = 'horizontal';
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Get the data
file_ID = fopen(whole_brain_file);
data = textscan(file_ID, '%f %f %f %f %f %f %f');
fclose(file_ID);
X = data{1}; % X coordinate in the altas
Y = data{2}; % Y coordiante in the atlas
Z = data{3}; % Z coordiante (Atlas slice)
region = data{6}; % The brain region the cell is assigned to.

% Exclude "cells" that are outside the brain
indexer = region==1;
disp(['Not including ' num2str(sum(indexer)) ' cells that are outside the brain.']);
region = region(~indexer);
X = X(~indexer);
Y = Y(~indexer);
Z = Z(~indexer);


% Update the region to the correct level and label all cells in regions
% with to few cells as 'other'.
[region, unique_regions, region_names] = HAN_update_level(region, level, min_cells);


% Load the brain template
disp('Loading the template volume...')
if strcmp(atlas_type,'50um atlas')
    load('Allen_3D_50um_template.mat','Allen_3D_template');
elseif strcmp(atlas_type, '10um atlas')
    load('Allen_3D_10um_template.mat','Allen_3D_template');
else
    error('Unknown atlas type')
end
disp('...done')
    

% Grab the colormap, add a line for the 'other' group
load('Allen_colormap.mat','cmap')
cmap = [cmap; [0.9 0.9 0.9]];


% Add the correct color to every cells
cell_color = zeros(length(region), 3);
for i=1:length(region)
    cell_color(i,:) = cmap(region(i), :);
end


% Now comes something fun, we're flipping X, Y and Z so we don't have to
% rewrite anything for different (e.g. sagittal) orientation
if strcmp(orientation, 'sagittal')
    temp = Z;
    Z = X;
    X = temp;
    
    Allen_3D_template = permute(Allen_3D_template, [1, 3, 2]);
end


% Same for horizontal slices
if strcmp(orientation, 'horizontal')
    temp = Z;
    Z = Y;
    Y = temp;
    
    Allen_3D_template = permute(Allen_3D_template, [3, 2, 1]);
end


% Find the unique atlas slices
atlas_slices = unique(Z);


% Deal with the fact that not all atlas slices are samples
[~, index] = min(abs(atlas_slices - start_slice));
start_slice = round(atlas_slices(index));


% Make the figure
main_figure = figure();


% Plot the template
template_projection = imagesc(squeeze(Allen_3D_template(:,:,start_slice)));
colormap gray
axis equal
hold on


% Plot the scatter on top
indexer = round(Z) == start_slice;
scatter(X(indexer), Y(indexer),1,cell_color(indexer,:));


flag = 0;
end

