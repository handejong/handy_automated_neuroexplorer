function [outputArg1,outputArg2] = HAN_find_cells(source_image_filename, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% Deal with input arguments
% Deal with input arguments
atlas_path = 'ARA_han/Slides ARA 10um/slide_';
atlas_type = '10um atlas';
save_data = false;
plot_data = false;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case 'plot data'
            plot_data = true;
            
        case 'save data' 
            save_data = true;
            
        case 'cell width'
            cell_width = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'min distance'
            min_distance = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end



outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

