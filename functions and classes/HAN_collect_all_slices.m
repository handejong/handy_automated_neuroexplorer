function whole_brain_data = HAN_collect_all_slices(output_filename, input_folder, include_slices, varargin)
%UNTITLED6 Summary of this function goes here
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
            
        case 'save data'
            save_data = true;
            
        case 'plot data'
            plot_data = true;
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Output data
whole_brain_data = [];

for slice = include_slices
    % Update the user
    disp(['Working on slice: ' num2str(slice)])
    
    % Path to the files
    cell_file = [input_folder num2str(slice) '_cells_mapped.txt'];
    
    % Get the cells
    file_ID = fopen(cell_file);
    cells_temp = textscan(file_ID, '%f %f %f %f %f %f');
    fclose(file_ID);
    
    % Convert the cells array to a normal matrix. In the last column put
    % the slice number we actually got these cells from.
    temp = [cells_temp{1} cells_temp{2} cells_temp{3} cells_temp{4} cells_temp{5} cells_temp{6} ones(size(cells_temp{6}))*slice];
    
    % Add to the final data
    whole_brain_data = [whole_brain_data; temp];
end


% If saving is on, save the data
if save_data
    dlmwrite(output_filename, whole_brain_data, 'delimiter', '\t');
end
  

% If plotting, plot data
if plot_data
    HAN_show_whole_brain(whole_brain_data,'level',level,'min cells',min_cells);
end



end

