function data = HAN_make_plots(whole_brain_file, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Deal with input arguments
plot_data = false;
save_data = false;
level = 1;
min_cells = 1;
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
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Grab the brain
% Get the data
file_ID = fopen(whole_brain_file);
temp_data = textscan(file_ID, '%f %f %f %f %f %f %f');
fclose(file_ID);
% X = temp_data{1};
% Y = temp_data{2};
% Z = temp_data{3};
region = temp_data{6};

% Exclude cells that are outside the brain
disp(['Not including ' num2str(sum(region==1)) ' cells that are outside the brain.']);
region = region(region~=1);


% Grab the colormap, add a line for the 'other' group
load('Allen_colormap.mat','cmap')
cmap = [cmap; [0.9 0.9 0.9]];


% Update the region to the correct level and label all cells in regions
% with to few cells as 'other'.
[region, unique_regions, region_names] = HAN_update_level(region, level, min_cells);


% Grab total cells
total_cells = length(region);

% Cells per region
for i=1:length(unique_regions)
    data(i).atlas_ID = unique_regions(i);
    data(i).region_name = region_names{i};
    data(i).cells = sum(region==unique_regions(i));
    data(i).cells_fraction = data(i).cells/total_cells;
end




% TODO, PUT CELL PER µm value

% Make X into a categorical
for i=1:length(data)
    X{i} = char(data(i).region_name);
    Y(i) = data(i).cells;
end

% A quirk of Matlab, but categoricals are alsays alphabetical, so you have
% to reorder them to get the original ordering.
X = reordercats(categorical(X),X);

% Plot the bars
main_figure = figure;
for i=1:length(data)
    bar(X(i),Y(i), 'FaceColor',cmap(data(i).atlas_ID,:));
    hold on
end
main_figure.Position(3) = 800;
ylabel('# of cells')
title('Total numer of cells')

end

