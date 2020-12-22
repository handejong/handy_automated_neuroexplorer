function flag = HAN_manual_inspector(source_image_filename, varargin)
%HAN_MANUAL_INSPECTOR Manually inspect and correct individual brain slices
%   Detailed explanation goes here
%
%   TODO

% Output flag
flag = 1;


% Some output text
disp(['Opening: ' source_image_filename])


% GLOBAL: Information about the mouse
c_point = []; % The current mouse point
mouser.empty = true; 
mouser.action = 'empty';
mouser.released = true;

% GLOBAL: Setting about what we are doing
new_cells_accepted = false;
currently_zoomed = false;
made_new_reference = false;

% GLOBAL: variables to make it possible to undo removal of cells
undo = []; undo{1} = [nan, nan]; % Removed cells in order to undo
undo_color = []; undo_color{1} = [nan nan nan]; % To sto
undo_all_data = []; undo_all_data{1} = [nan nan nan nan nan nan];

% Other Global data
output_reference = [];
training_data = []

% Deal with input arguments
edited = false;
save_data = false;
static = false;
mapped = false;
level = 3;
min_cells = 1;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case 'edited'
            edited = true;
            
        case 'mapped'
            mapped = true;
            
        case 'level'
            level = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'min cells'
            min_cells = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'save data' 
            save_data = true;
            
        case 'static' % only visualization, no edditing
            static = true;
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Grab all the paths
paths = HAN_get_paths(source_image_filename);


% Import the source image
source_image = imread(source_image_filename);
source_dims = size(source_image);


% figure out the image type
switch class(source_image)
    
    case 'uint8'
        image_type = 2^8;
        
    case 'uint16'
        image_type = 2^16;
        
    case 'uint32'
        image_type = 2^32;
    
    otherwise
        disp(['Unknown class ''' class(source_image) '''but assuming 16bit'])
        max_image_brightness 
end


% If the source image is grayscale, make it RGB anyway, so we won't have
% trouble plotting the atlas in color on top
if length(source_dims)<3
    source_image = repmat(source_image,[1 1 3]);
end


% Another manipulation of the source image is that it's common to take
% microscopy images somewhat at the low end of the dynamic range, so they
% might appear dark on the screen. If the picture is really dark, we'll
% increase the brightness somewhat.
temp_max_value = max(source_image(:));
if temp_max_value<image_type
    source_image = source_image.* (image_type / temp_max_value);
end


% Get the cells
if edited
    disp('Loading the manually edited cell file.')
    file_ID = fopen(paths.cells_edited);
    cells_temp = textscan(file_ID, '%f %f ', 'HeaderLines', 3);
    fclose(file_ID);
    cells(:,1) = cells_temp{1};
    cells(:,2) = cells_temp{2};
elseif mapped
    disp('Loading cells that are mapped to the reference atlas.')
    file_ID = fopen(paths.cells_mapped);
    cells_temp = textscan(file_ID, '%f %f %f %f %f %f');
    fclose(file_ID);
    cells(:,1) = cells_temp{4};
    cells(:,2) = cells_temp{5};
    region = cells_temp{6};
    all_data = [cells_temp{1}, cells_temp{2}, cells_temp{3}, cells_temp{4}, cells_temp{5}, cells_temp{6}];
else
    disp('Loading cells.')
    file_ID = fopen(paths.cells);
    cells_temp = textscan(file_ID, '%f %f', 'HeaderLines', 3);
    fclose(file_ID);
    cells(:,1) = cells_temp{1};
    cells(:,2) = cells_temp{2};
end

% If mapped cells, we also want to know the cell color
if mapped
    
    % Update the region to the correct level and lavel cells in regions
    % with to few cells as 'other'.
    [region, unique_regions, region_names] = HAN_update_level(region, level, min_cells);
    
    % Grab the colormap, add a line for the 'other' group
    load('Allen_colormap.mat','cmap')
    cmap = [cmap; [0.9 0.9 0.9]];

    % Add the correct color to every cell
    cell_color = zeros(length(region), 3);
    for i=1:length(region)
        cell_color(i,:) = cmap(region(i), :);
    end
end


% Make a backup
cells_original = cells;

% Get the reference
reference = imread(paths.ref);

% Get the AIDAHist_setings
settings_updated = false;
load(paths.AIDAHis_settings, 'settings')
settings; %<- just remember that settings means AIDAHist_settings

    
% Make the main figure
main_figure = figure('Position',[490, 200, 1000, 800]);


% Plot the source image
source_projection = image(source_image,...
    'Tag','source'); hold on;


% Plot the cells
if mapped
    % plot the color as well
    cell_scatter = scatter(cells(:,1), cells(:,2), 1, cell_color,...
        'Tag','cell_scatter');
else
    cell_scatter = scatter(cells(:,1), cells(:,2), 'x','r',...
        'Tag','cell_scatter');
end
axis equal
original_dims = [xlim(); ylim()];


% Plot the reference, make it invisible
reference_projection = imagesc(reference,...
    'Visible','off',...
    'AlphaData',ones(size(reference)).*0.2,...
    'Tag','ref');
caxis([0 1000])


% Add callbacks
if ~static
    main_figure.KeyPressFcn = @key_press;
    reference_projection.ButtonDownFcn = @mouse_down;
    source_projection.ButtonDownFcn = @mouse_down;
    cell_scatter.ButtonDownFcn = @mouse_down;
    main_figure.WindowButtonMotionFcn = @mouse_move;
    main_figure.WindowButtonUpFcn=@mouse_release;
end
set(gcf,'Visible','on')


% Make the function run untill the user is finished
we_are_done = false;
while ~we_are_done
    % Keep running untill the main figure is closed
    if ~ishandle(main_figure)
        we_are_done = true;
    end
    
    pause(0.01);
end


% If the user has not accepted the new cells, make sure this is on purpose
if ~new_cells_accepted && save_data
     input = questdlg('Accept the new cells?','Cell selection not accepted.','yes','no','no');
     if strcmp(input, 'yes')
         new_cells_accepted = true;
     else
         new_cells_accepted = false; %still
     end
end


% Update the cells if we do
if new_cells_accepted && save_data
    fileID = fopen(paths.cells_edited,'w');
    fprintf(fileID, 'AIDAhisto - Number of detected cells (EDDITED): %i\ncell postions (xy):\n\n',length(cells));
    for i=1:length(cells)
        fprintf(fileID,'%i %i\n', cells(i,1), cells(i,2));
    end
    fclose(fileID);
    disp(['Saved cells to: ' paths.cells_edited])
    
    % also updating the cells_mapped file if we are using that one
    if mapped
        % get the file path
        filename = paths.cells_mapped; 
        dlmwrite(filename, all_data, 'delimiter', '\t');
        disp(['Saved mapped cells to: ' paths.cells_mapped])
        
        % We might save a new reference as well
        if made_new_reference
            disp('Saving updated reference and training data.')
            imwrite(output_reference, paths.ref);
            save(paths.training_data, 'training_data');
            HAN_update_atlas_data(paths.source, training_data)
        end
    end
    
    % also save the AIDAHis_settings
    if settings_updated
        disp(['Saving new AIDAhis settings to: ' paths.AIDAHis_settings]);
        save(paths.AIDAHis_settings, 'settings');
    end
    
else
    warning('nothing saved')
end


% OK we made it
flag = 0;
disp('Done')


%%%%%%%%%%%%%%%%%%%%%% sub functions and callbacks %%%%%%%%%%%%%%%%%%%%%%%%
    
    function key_press(scr, ev)
    % Controls responses to key presses

    switch ev.Key


        % Manipulating the source image
        case 'b' %Increase the brightness
            source_projection.CData = source_projection.CData + 0.05*image_type;

        case 'v' %Decrease the brightness
            source_projection.CData = source_projection.CData - 0.05*image_type;

        case 'c' %Increase contrast
            source_projection.CData = source_projection.CData.*2;

        case 'x' %Decrease contrast
            source_projection.CData = source_projection.CData./2;

        case 'l' %Flip left-right
            flip_source_horizontal();

        case 'o' %Reset original
            source_projection.CData = source_image;
            
        
        % Manipulating the reference
        case 's' % Show the reference
            if strcmp(reference_projection.Visible, 'on')
                reference_projection.Visible = 'off';
            else
                reference_projection.Visible = 'on';
            end
            
            
        % Manipulating the cells
        case 'space'% Show the identified cells
            if strcmp(cell_scatter.Visible, 'on')
                cell_scatter.Visible = 'off';
            else
                cell_scatter.Visible = 'on';
            end
            
        case 'u' %undo last removal
            if length(undo) > 1
                cells = [cells; undo{end}];
                undo = undo(1:end-1);
                cell_scatter.XData = cells(:,1);
                cell_scatter.YData = cells(:,2);
                
                if mapped %more to undo
                    cell_color = [cell_color; undo_color{end}];
                    cell_scatter.CData = cell_color;
                    all_data = [all_data; undo_all_data{end}];
                    undo_color = undo_color(1:end-1);
                    undo_all_data = undo_all_data(1:end-1); 
                end
                
            end
            
        case 'r' %revert to original
            cells = cells_original;
            cell_scatter.XData = cells(:,1);
            cell_scatter.YData = cells(:,2);
            
        case 'return' % Accept these cells and close
            new_cells_accepted = true;
            delete(main_figure)
            we_are_done = true;
            disp('New cells accepted.')
        
            
        % Changing the threshold and try to find new cells
        case 'p'
            if ~mapped
                settings_updated = true;
                input = inputdlg("New Threshold:", "Give the new threshold.");
                undo{end+1} = cells;
                settings.thresw = str2double(input{1});
                [~, cells, ~] = AIDAhisto(paths.source,...
                    'SETTINGS',settings);
                cell_scatter.XData = cells(:,1);
                cell_scatter.YData = cells(:,2);
            else
                warning('Redoing the thresholding is currently not possible when inspecting mapped cels');
            end
         
            
        % Updating the reference
        case 'q'
            update_reference();
        
            
        % Unboud key
        otherwise
            disp(['Unbound key: ' ev.Key])
    end


end

    function mouse_down(scr, ev)
       % Deals with mouse clicks
       
       % The mouse button is still clicked
       mouser.released = false;
       
       % RIGHT mouse click -> standard action (remove cells)
       if ev.Button==3
           mouser.action = 'standard';
           temp_color = [0.8 0.4 0];
       end
        
       % LEFT mouse click -> zoom mode
       if ev.Button==1
           if ~currently_zoomed
               mouser.action = 'zoom';
               temp_color = [1 1 1];
           else % zoom out
               currently_zoomed = false;
               xlim(original_dims(1,:));
               ylim(original_dims(2,:));
               mouser.action = 'empty';
               mouser_release = true;
               return
           end
       end
       
       % Figure out what was actually clicked upon
       switch scr.Tag % what was clicked?
                
            case {'ref', 'source', 'cell_scatter'} % on the image
                if ~mouser.empty % there is no rectangle currently
                    delete(mouser.drawing)

                end

                mouser.start = [c_point(1), c_point(2)];
                mouser.drawing = rectangle(...
                    'Position',[c_point(1), c_point(2), 1, 1],...
                    'EdgeColor', temp_color);

           otherwise
                mouser.action = 'empty';
                mouser.released = true;
                disp('no idea what the mouse is doing')
        end

    end

    function mouse_move(scr, ev)
        % Deals with mouse movement
        
        
        % Get the current point
        punter=get(scr.Children,'CurrentPoint');
        c_point = punter(1,1:2);
        
        % Only do the rest if the mouse is not release
        if mouser.released; return; end
        
        
        % Otherwise do things
        switch mouser.action
            
            % Manipulate the rectangle
            case {'standard', 'zoom'}
                block=[mouser.start(1), mouser.start(2), c_point(1)-mouser.start(1), c_point(2)-mouser.start(2)];

            otherwise
                X = 5; % Mouse is not doing any action
                return
        end
        
        % Deal with rectangles drawn up or to the left.
        if block(3)<=0
            block(1)=block(1)+block(3);
            block(3)=block(3)*-1;   
        end
        if block(4)<=0
            block(2)=block(2)+block(4);
            block(4)=block(4)*-1;
        end
        mouser.drawing.Position = block;
        

    end

    function mouse_release(scr, ev)
        % Deals with mouse motion end
        
        mouser.released = true;
        
        switch mouser.action
            
            case 'standard'
                selection = round(mouser.drawing.Position); %[X Y W H]
                delete(mouser.drawing)
                mouser.action = 'empty';
                
                % Delete all the data points in the rectangle
                X = [selection(1), selection(1), selection(1)+selection(3), selection(1)+selection(3)];
                Y = [selection(2), selection(2)+selection(4), selection(2)+selection(4), selection(2)];
                indexer = inpolygon(cells(:,1), cells(:,2), X, Y);
                
                % Allow for this to be undone by saving the to-be-removed
                % cells
                undo{end+1} = cells(indexer, :);
                
                % If working with mapped data, there is more to undo
                if mapped
                    undo_color{end+1} = cell_color(indexer, :);
                    undo_all_data{end+1} = all_data(indexer, :);
                end
                
                % remove them from the cells
                cells = cells(~indexer, :);
                
                % Update the plot
                cell_scatter.XData = cells(:,1);
                cell_scatter.YData = cells(:,2);
                
                % If working on mapped cells, update a bit more
                if mapped
                    cell_color = cell_color(~indexer, :);
                    cell_scatter.CData = cell_color;
                    all_data = all_data(~indexer, :);
                end
                
            case 'zoom'
                selection = round(mouser.drawing.Position); %[X Y W H]
                delete(mouser.drawing)
                mouser.action = 'empty';
                
                % Zoom to the rectangle
                xlim([selection(1), selection(1) + selection(3)])
                ylim([selection(2), selection(2) + selection(4)])
                
                % Currently zoomed
                currently_zoomed = true;
                
            otherwise
                X= 5; %do nothing
        end

    end
    
    function flip_source_horizontal()
        %UNDOCUMENTED function that will flip the source image INCLUDING
        %THE IMAGE FILE ITSELF horizontally. Should actually be done during
        %preprocessing of the dataset, this function is for emergencies
        %only.
        
        disp('hoi')
        
        % Flip the image
        source_image = fliplr(source_image);
        source_projection.CData = source_image;
        imwrite(source_image, paths.source);        
        
        % Flip the reference
        reference = fliplr(reference);
        reference_projection.CData = reference;
        imwrite(reference, paths.ref);
        
        % Flip the cells (note that the original file is not changed)
        cells(:,1) = -1*(cells(:,1) - 0.5 * size(source_image, 2)) + 0.5*size(source_image, 2);
        cell_scatter.XData = cells(:,1);
        cell_scatter.YData = cells(:,2);
        
    end

    function update_reference()
        %UPDATE_REFERENCE opens the make reference GUI and let's the user
        %map a new reference to the selected image
        
        
        % Figure out the started slice
        start_slice = 200;
        if mapped; start_slice = cells_temp{3}(1); end
        
        % Figure out what atlas to map to
        if isfile(paths.atlas_data)
            load(paths.atlas_data,'atlas_data');
            atlas_type = atlas_data(1).atlas_type;
        else
            % Ask the user
            atlas_type = questdlg('Map to which atlas?','Atlas type',...
                '50um atlas', '10um atlas', '10um atlas'); 
        end
        
        % Make the reference
        [output_reference, training_data] = HAN_make_reference(paths.source,...
            'start slice',start_slice,...
            atlas_type);
        
        % Check if the user accepted the map, otherwise just return
        if output_reference ==0; return; end
        
        % Update that the cells are now mapped
        mapped = true;
        made_new_reference = true;
        
        % Update the reference
        reference_projection.CData = output_reference;
        
        % If save data, mapp the cells to the reference
        all_data = HAN_map_cells_to_atlas(paths.source,...
            'use cells', cells,...
            'use training data', training_data);
        
        % Update the cell color and all data
        [region, unique_regions, region_names] = HAN_update_level(all_data(:,6), level, min_cells);

        % Add the correct color to every cell
        cell_color = zeros(length(region), 3);
        for i=1:length(region)
            cell_color(i,:) = cmap(region(i), :);
        end
        cell_scatter.CData = cell_color;
        
    end

end

