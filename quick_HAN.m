%% Start
close all
clear all
clc
addpath('AIDAhisto_Matlab') % Rewriten version of AIDAhisto (Original code by Niklas Pallast).
addpath('functions and classes') % All other code the neuroexplorer uses
addpath('Allen Reference Atlas') % The modified reference atlas and its legend


%% Variables
path_to_the_brain = '/Users/handejong/Box Sync/Temp/H108_0 brain slices jpg/';
file_ext = '.jpg';
atlas_type = '50um atlas';
automatic_reference = false; % Should we automatically allign the atlas or open a GUI to do it manually
slice_numbers = [10:20];
atlas_thickness = 50; % thickness of the atlas slices in µm
cell_threshold_value = 1.5; % The threshold above which cells are identified is automatically calculated, but you can 'weigh' it here. Use trial and error.
cell_size = 12; % How many pixels a cell is wide (on average) in your pictures (should be at least 3).
channel_number = 1; % If your source images are multi-chanel images (for instance RGB = red, green, blue) select the channel you want to do the analysis in here.
reference_channel = 2; % You can put 0 here if you don't want to use an artifact channel, or only have grayscale images.
save_data = true; % Save the data or just quickly visualize
plot_data = true; % Show data after each analysis
start_slice = 201; % The first reference atlas slice that will be shown
slice_thickness = 100; % In µm, used to predict the next most likely atlas slice and to figure out if any slices are missing
new_filename = 'Han test.txt';
level = 4;
min_cells = 100;


%% Go through all the slices
for slice = slice_numbers
    
     % Grab the filename
     full_file = [path_to_the_brain num2str(slice) file_ext];
     
     % Output text
     disp(['Working on slice: ' full_file]);
     
     % Figure out if we should save the data
     save_str = 'skipp'; if save_data; save_str='save data'; end
     plot_str = 'skipp'; if plot_data; plot_str='plot data'; end
     
     % AUTOMATIC reference
     if automatic_reference
         [reference, training_data] = HAN_automatic_reference(full_file, 'start slice', start_slice',...
             save_str, plot_str, '50um atlas');
     else
         [reference, training_data] = HAN_make_reference(full_file,save_str,atlas_type);
     end
     
     % Find cells
     fprintf("Processing wholebrain slice...\n")
     AIDAhisto(full_file,...
        'WIDTH', cell_size,...
        'DARK_PEAKS', 1,...
        'BAR_FILTER', 1,...
        'CHANNEL',channel_number,...
        'REF_CHANNEL',reference_channel,...
        'THRES_W', cell_threshold_value,...
        'SAVE_DATA', save_data,...
        'PLOT_DATA', plot_data);
    
    % Map to refference
    HAN_map_cells_to_atlas(full_file, save_str, plot_str);
    
    % Grab the atlas slice for the next run
    start_slice = training_data.ref_nr - round(slice_thickness/atlas_thickness);

end


%% Run the manual inspecter on every slide
for slice = slice_numbers
    % Output text
    disp(['Working on slice ' num2str(slice)])
    
    % Grab the filename
    full_file = [path_to_the_brain num2str(slice) file_ext];
    
    % Open the manual inspector
    flag = HAN_manual_inspector(full_file, 'save data', 'mapped');
end


%% Collect all the data in one file


% Figure out if we should plot and/or save data
save_str = 'skipp'; if save_data; save_str='save data'; end
plot_str = 'skipp'; if plot_data; plot_str='plot data'; end

% Output data
HAN_collect_all_slices(new_filename, path_to_the_brain, slice_numbers, save_str, plot_str, 'level', level);


%% Show the plots and the 3D brain
% Show the entire brain in 3D
HAN_show_whole_brain(new_filename,...
    'level', level,...
    'min cells', min_cells);

% Show some plots
data = HAN_make_plots(whole_brain_file, 'level', 5, 'min cells', 100);
    