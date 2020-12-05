function [updated_region, unique_regions, region_names] = HAN_update_level(region, level, min_cells)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% Note that region 1328 means 'other'


% Find unique brain regions
unique_regions = unique(region);

% Find the brain region name as well as the atlas ID at a certain level
disp('Looking up brain region names....')
[~, new_id] = HAN_look_up_Allen(unique_regions, level);

% Relabel all the cells to the region_id at that level
for i=1:length(unique_regions)
    index = region == unique_regions(i);
    region(index) = new_id(i);
end

% Plot only the regions with a minimum number of cells in there
unique_regions = unique(region);
for i=1:length(unique_regions)
    index = region == unique_regions(i);
    if sum(index(:))<min_cells
        new_id(i) = 1328; % means 'other'
        region(index) = new_id(i);
    end
end

% Ok and find the unique regions again
unique_regions = unique(region);
region_names = HAN_look_up_Allen(unique_regions, level);
disp('... done')


% Output arguments
updated_region = region;

end

