function flag = HAN_find_missing_slides(whole_brain_file, slice_thickness, atlas_slice_thickness, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Start flag
flag = 1;

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

% Get the data
file_ID = fopen(whole_brain_file);
data = textscan(file_ID, '%f %f %f %f %f %f %f');
fclose(file_ID);
X = data{1};
Y = data{2};
Z = data{3};
microscopy_slides = data{7};
region = data{6};

% Unique slices
unique_slices = unique(microscopy_slides);
slice_to_atlas = [unique_slices, zeros(size(unique_slices))];

% Figure out which atlas slice was mapped to which slice
for i=1:length(slice_to_atlas)
    index = find(slice_to_atlas(i,1)==microscopy_slides);
    slice_to_atlas(i,2) = Z(index(1));
end

% Overview of all the atlas slices
[atlas_slices, ind] = unique(slice_to_atlas(:,2));
duplicate_ind = setdiff(1:length(slice_to_atlas(:,2)), ind);

% Print out overview of double sampled multiple times
disp(' ');
disp('The following slices are sampled double, you might want to redo those:')
for i=1:length(duplicate_ind)
    a_slice = slice_to_atlas(duplicate_ind(i),2);
    s_slice = slice_to_atlas(duplicate_ind(i),1);
    temp = find(slice_to_atlas(:,2)==a_slice);
    s1_slice = slice_to_atlas(temp(1),1);
    disp(['Slice ' num2str(s1_slice) ...
        ' and slice ' num2str(s_slice) ...
        ' are both mapped to: ' num2str(a_slice) '.']);
end
if isempty(duplicate_ind); disp('none'); end

% Print how many slices might be missing
total_a_distance = (max(atlas_slices) - min(atlas_slices)) * atlas_slice_thickness + atlas_slice_thickness;
total_s_distance = length(slice_to_atlas(:,1)) * slice_thickness;
missing_distance = total_a_distance - total_s_distance;
missing_slices = round(missing_distance/slice_thickness);
disp(' ')
if missing_distance>0
    disp(['Missing data for ' num2str(missing_distance) 'µm (' num2str(missing_slices) ' slices).']);
    disp(['Microscopy distance covered: ' num2str(total_s_distance) 'µm'])
    disp(['Atlas distance covered: ' num2str(total_a_distance) 'µm'])
else
    disp(['Probably some misallignment or variation in slice thickness:'])
    disp(['Microscopy distance covered: ' num2str(total_s_distance) 'µm'])
    disp(['Atlas distance covered: ' num2str(total_a_distance) 'µm'])
end




% Finished
flag = 0;
end

