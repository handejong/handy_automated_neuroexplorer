function [output_reference, training_data] = HAN_make_reference(source_image_filename, varargin)
%HAN_MAKE_REFERENCE produces ARA brain atlas reference images to microscopy
%images (input) of the mouse brain.

%   The input image is a microscopy image of a coronal slice of a mouse
%   brain. 'Output_reference' is the output reference atlas picture
%   transformed form the slides in the ARA. Training data is a .mat file
%   containing a converted input image, the output reference and the
%   landmarks. The training_data can be used to train a neural network to
%   automate this proces.
%
%   HAN_make_reference is mostly controlled by the keyboards. The possible
%   keys are printed to the command line at startup.
%
%   HAN_make_reference is part of the Handy Automated Neuroexplorer (HAN).
%   HAN is made by Johannes de Jong, j.w.dejong@berkeley.edu.


% Global variables
c_atlas_slice = 200; % Current atlas slice
landmark_mode = false; % Landmark mode on or off
show_output_reference = false; % Showing the output reference
showing_landmark_numbers = false; % Showing landmark numbers
atlas_dims = []; % Atlas dimensions
source_dims = []; % Source dimensions
landmark_number_handles = {};
tform =[]; % Transformation affine2d object


% Mouse and keyboard action
c_point = [];
mouse_target = [];


% X and Y coordinates of landmarks
atlas_landmarks = [nan nan];
source_landmarks = [nan nan];


% Startup text
disp(' ')
disp('...')
disp('To use the GUI, use the following keypresses.')
disp('Arrows: to select the appropriate atlas slice. ')
disp('Space: to go into landmark mode.')
disp('d or f: to delete 1 or all landmarks.')
disp('n: to show the landmark numbers.')
disp('t: to perform the atlas transformation')
disp('s: to show the transformation.')
disp(' ')
disp('To adjust the contrast and brighness of the source image use:')
disp('c and x: to increase/decrease contrast')
disp('b and v: to increase/decrease brightness')
disp('o: to reset tue original source image')
disp(' ')
disp('Finally, use ''return'' to accepted the transformed reference and close the GUI to finish.')
disp(' ')


% Deal with input arguments
atlas_path = 'Allen Reference Atlas/Slides 10um/slide_';
atlas_type = '10um atlas';
save_data = false;
for i=1:length(varargin)
    switch varargin{i}
        
        case 'skipp'
            continue;
            
        case '50um atlas' % Select the 50um atlas
            atlas_path = 'Allen Reference Atlas/Slides 50um/slide_';
            atlas_type = '50um atlas';
            
        case '10um atlas' % Select the 10um atlas
            % This is the default, so doesn't change anyting for now
            atlas_path = 'Allen Reference Atlas/Slides 10um/slide_';
            atlas_type = '10um atlas';
            
        case 'start slice' % set the atlas start slice
            c_atlas_slice = varargin{i+1};
            varargin{i+1} = 'skipp';
            
        case 'save data' 
            save_data = true;
            
        otherwise
            disp('Unknown input argument:')
            varargin{i};
    end
end


% Error handeling 1, can we find the atlas?
if ~isfile([atlas_path '1.tiff'])
    error('Unable to locate the Allen reference atlas.')
end


% Figure out the highest numbered atlas slide
if strcmp(atlas_type, '10um atlas')
    max_atlas_slice = 1320; % Highest numbered slice in the atlas
    slide_used = false(1320,1);
elseif strcmp(atlas_type, '50um atlas')
    max_atlas_slice = 264;
    slide_used = false(264,1);
end


% Error handeling 2, if the user wants a slide that is not in the atlas
if c_atlas_slice>max_atlas_slice || c_atlas_slice < 1
    warning(['Slide ' num2str(c_atlas_slice) ' is not in the atlas.'])
    c_atlas_slice = 200;
end


% Grab all the paths
paths = HAN_get_paths(source_image_filename);

% See if there is allready atlas data
if isfile(paths.atlas_data)
    disp('Loading atlas data.')
    load(paths.atlas_data, 'atlas_data')
else
    atlas_data = [];
end

% Check form the atlas data if any slide is allready used
for i=1:length(atlas_data)
    slide_used(atlas_data(i).ref_nr) = true;
end

% Load the input image
source_image = imread(paths.source);
source_dims = size(source_image);
% ... error handeling


% Load the first reference atlas image
atlas_image = imread([atlas_path num2str(c_atlas_slice) '.tiff']);
atlas_dims = size(atlas_dims);


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

% Grab the colormap, add a line for the 'other' group
load('Allen_colormap.mat','cmap')
cmap = [cmap; [0.9 0.9 0.9]];

% Output arguments
output_reference = zeros(source_dims(1:2));


% Make the main figure
main_figure = figure('Position', [1,1, 1600, 800],...
    'Units','Normalized');


% Plot the atlas
plot_1 = subplot(1,2,1);
atlas_projection = imagesc(atlas_image,...
    'Tag','atlas');
axis equal; axis off;
colormap(cmap);


% Show a little warning if the slide is allready used
slide_used_warning = text(10, 10, 'Slide allready mapped',...
    'Visible','off');


% Update the warning if the slide is allready mapped
if slide_used(c_atlas_slice)
    slide_used_warning.Visible = 'on';
else
    slide_used_warning.Visible = 'off';
end


% Print the slide number in the other corner
slide_number_text = text(210, 10, num2str(c_atlas_slice));


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


% Plot the source image (always in color)
plot_2 = subplot(1,2,2);
image_projection = image(source_image,'Tag','source');
axis equal; axis off; hold on;


% Plot the output reference on top
output_reference_projection = imagesc(output_reference,...
    'Visible','off',...
    'AlphaData',ones(size(output_reference)).*0.2,...
    'Tag','ref');
colormap(cmap);


% Work on the the GUI
main_figure.Children(2).Position = [0.0, 0.05, 0.5, 1];
main_figure.Children(1).Position = [0.5, 0.05, 0.5, 1];


% Plot the landmarks on top, but make them invisible
subplot(plot_1); hold on;
atlas_landmark_plot = scatter([], [],...
    'Visible','off');
subplot(plot_2); hold on;
source_landmark_plot = scatter([], [],...
    'Visible','off');


% Add callbacks 
main_figure.KeyPressFcn = @key_press;
atlas_projection.ButtonDownFcn = @mouse_down;
image_projection.ButtonDownFcn = @mouse_down;
output_reference_projection.ButtonDownFcn = @mouse_down;
main_figure.WindowButtonMotionFcn = @mouse_move;
set(gcf,'Visible','on')


% Make the function run untill the user is finished
we_are_done = false;
reference_accepted = false;
while ~we_are_done
    % Keep running untill the main figure is closed
    if ~ishandle(main_figure)
        we_are_done = true;
    end
    
    pause(0.01);
end


% If the user has not accepted the reference, make sure this is on purpose
if ~reference_accepted
     input = questdlg('Accept this reference?','Reference not accepted.','yes','no','no');
     if strcmp(input, 'yes')
         reference_accepted = true;
     else
         reference_accepted = false; %still
     end
end
    

% Finish up
if ~reference_accepted
    output_reference = 0;
    training_data = 0;
    return
else
    training_data.source = source_image;
    training_data.atlas = atlas_type;
    training_data.ref_nr = c_atlas_slice;
    training_data.atlas_landmarks = atlas_landmarks;
    training_data.source_landmarks = source_landmarks;
    training_data.tform = tform;
end


% Save the data if requested
if save_data
    imwrite(output_reference, paths.ref);
    save(paths.training_data, 'training_data');
    
    % Update the atlas data
    HAN_update_atlas_data(source_image_filename, training_data);
end


%%%%%%%%%%%%%%%%%%%%%%% sub functions and callbacks %%%%%%%%%%%%%%%%%%%%%%%


    function key_press(scr, ev)
        % Controls responses to key presses
        
        switch ev.Key
            
            % Navigating the atlas
            case 'rightarrow'
                c_atlas_slice = max([c_atlas_slice-1, 1]);
                update_atlas();
                
            case 'leftarrow'
                c_atlas_slice = min([c_atlas_slice+1, max_atlas_slice]);
                update_atlas();
                
            case 'downarrow'
                c_atlas_slice = min([c_atlas_slice+25, max_atlas_slice]);
                update_atlas();
                
            case 'uparrow'
                c_atlas_slice = max([c_atlas_slice-25, 1]);
                update_atlas();
            
                
            % Manipulating the source image
            case 'b' %Increase the brightness
                image_projection.CData = image_projection.CData + 0.05*image_type;
                
            case 'v' %Decrease the brightness
                image_projection.CData = image_projection.CData - 0.05*image_type;
                
            case 'c' %Increase contrast
                image_projection.CData = image_projection.CData.*2;
                
            case 'x' %Decrease contrast
                image_projection.CData = image_projection.CData./2;
                
            case 'l' %Flip left-right
                flip_source_horizontal();

            case 'o' %Reset original
                image_projection.CData = source_image;
            
                
            % About the landmarks    
            case 'space' %Landmark mode
                landmark_mode = ~landmark_mode;
                update_landmarks();
                
            case 'd' %Delete last landmark
                source_landmarks = source_landmarks(1:end-1,:);
                atlas_landmarks = atlas_landmarks(1:end-1,:);
                update_landmarks();
                
            case 'f' % Start over
                source_landmarks = [nan nan];
                atlas_landmarks = [nan nan];
                update_landmarks();
                
            case 'n' %Show landmark numbers
                show_landmark_numbers();
                
            case 'a' %Auto landmarks
                auto_landmarks();
            
                
            % Work out the output reference
            case 't' % Transform the selected reference
                show_output_reference = true;
                output_reference = transform_reference();
                update_atlas();
                
            case 's' %Show the output reference
                show_output_reference = ~show_output_reference;
                update_atlas();
                
            case 'return' %Reference accepted
                reference_accepted = true;
                disp('Reference accepted!')
                delete(main_figure);
                
            % Unboud key
            otherwise
                disp(['Unbound key: ' ev.Key])
        end
        
        
    end

    function mouse_down(scr, ev)
        % Deals with mouse clicks
        
        % If landmark mode, make a new landmark
        if landmark_mode
            
            % Only add landmark if mouse is on target
            if mouse_on_target
                if strcmp(mouse_target, 'atlas')
                   atlas_landmarks = [atlas_landmarks; c_point];
                else
                   source_landmarks = [source_landmarks; c_point];
                end
                
                % And update
                update_landmarks();
            end
            
        end

    end

    function mouse_move(scr, ev)
        % Deals with mouse movement
        
        % Get the current point
        punter=get(scr,'CurrentPoint');
        %c_point = punter(1,1:2);

        % Figure out what was clicked on and get the relative c_point;
        if punter(1)>0.5
            mouse_target = 'source';
            punter=get(scr.Children(1),'CurrentPoint');
            c_point = punter(1,1:2);
            
        else
            mouse_target = 'atlas';
            punter=get(scr.Children(2),'CurrentPoint');
            c_point = punter(1,1:2);
            
        end
        
        
    end

    function update_atlas()
        % Updates the atlas

        % Update the atlas slice
        filename = [atlas_path num2str(c_atlas_slice) '.tiff'];
        atlas_projection.CData = imread(filename);
        
        % Update the transformed (warped) atlas
        if show_output_reference
            output_reference_projection.Visible = 'on';
            output_reference_projection.CData = output_reference;
            output_reference_projection.AlphaData = ones(size(output_reference)).*0.2;
            output_reference_projection.AlphaData(output_reference==1)=0;
        else
            output_reference_projection.Visible = 'off';
        end
        
        % Update the warning if the slide is allready mapped
        if slide_used(c_atlas_slice)
            slide_used_warning.Visible = 'on';
        else
            slide_used_warning.Visible = 'off';
        end
        
        % Update the atlas slice number
        slide_number_text.String = num2str(c_atlas_slice);
        
        % draw
        drawnow(); 
        
    end

    function update_landmarks()
        % Updates the landmarks
        
        % Update the landmark plots
        atlas_landmark_plot.XData = atlas_landmarks(:,1);
        atlas_landmark_plot.YData = atlas_landmarks(:,2);
        source_landmark_plot.XData = source_landmarks(:,1);
        source_landmark_plot.YData = source_landmarks(:,2);
        
        % Remove nan landmarks (only after the first landmark was created
        if size(atlas_landmarks, 1) == 2
            atlas_landmarks = atlas_landmarks(~isnan(atlas_landmarks(:,1)),:);
        end
        if size(source_landmarks, 1) == 2
            source_landmarks = source_landmarks(~isnan(source_landmarks(:,1)),:);
        end
        
        % Plot or do not plot the landmarks
        if landmark_mode
            atlas_landmark_plot.Visible = 'on';
            source_landmark_plot.Visible = 'on';
        else
            atlas_landmark_plot.Visible = 'off';
            source_landmark_plot.Visible = 'off';
        end
        
        % draw
        drawnow(); 

    end

    function bool = mouse_on_target()
        % Checks if the mouse is really on the target
        
        
        % Smaller then 0 is of target anyway, let's speed things up
        if min(c_point)<=0
            bool = false;
            return
        end
        
        % Find the target size
        if strcmp(mouse_target, 'atlas')
            max_x = atlas_dims(2);
            max_y = atlas_dims(1);
        else
            max_x = source_dims(2);
            max_y = source_dims(1);
        end
        
        % Bigger then the max x or max y -> also off target
        if c_point(1)>max_x || c_point(2)>max_y
            bool = false;
        end
        
        % Apparently on target
        bool = true;
        
    end
    
    function reference = transform_reference
        % Transforms the reference based on the landmarks
        
        % Check if the number of landmarks on the atlas and the source
        % picture actually correspond.
        if length(atlas_landmarks) ~= length(source_landmarks)
            errordlg('The atlas and the source image do not have the same number of landmarks. Use ''d'' to delete landmarks or ''f'' to start over. Use ''n'' to show the landmark numbers.')
            reference = output_reference;
            return;
        end
        
        % Use the landmarks to find the tform matrix
        [tform,inlierPtsDistorted,inlierPtsOriginal] = ...
        estimateGeometricTransform(atlas_landmarks, source_landmarks,...
        'affine','MaxDistance',200);
    
        % Include some presentation of the points that were actually used
        % ....
 
        % I'd like to point out that the transformation between two sets of
        % landmarks is simply done by solving for theta H1 = H0*theta + b;
        % The easiest way to do this is simply to add a column of 1s to the
        % landmark matrixes and then to solve using:
        %
        %   >>t_form = atlas_landmarks\source_landmarks;
        %
        % ALL THIS OTHER CODE IS JUST FOR THE VISUALIZATION YOU DON'T NEED
        % IT. Altough it is nice that estimateGeometricTransform trows out
        % outliers caused by strongly deformed microscopy slices or simply
        % the user now knowing their anatomy.
        
        % Grab the target image
        filename = [atlas_path num2str(c_atlas_slice) '.tiff'];
        target_image = imread(filename);
        
        % Do the transformation on the image
        outputView = imref2d(size(source_image));
        reference = imwarp(target_image,tform,'OutputView',outputView);
    end
    
    function show_landmark_numbers()
        % Shows landmark numbers, or removes them

        if ~showing_landmark_numbers
            % Empty out old landmark numbers
            landmark_number_handles = {};
            counter = 1;
            
            subplot(plot_1)
            for j=1:size(atlas_landmarks, 1)
                landmark_number_handles{counter} = text(atlas_landmarks(j,1)+2,...
                    atlas_landmarks(j,2), num2str(j));
                counter = counter+1;
            end
                
            subplot(plot_2)
            for j=1:size(source_landmarks, 1)
                landmark_number_handles{counter} = text(source_landmarks(j,1)+2,...
                    source_landmarks(j,2), num2str(j));
                counter = counter+1;
            end
            
            % Mark landmark numbers as being shown.
            showing_landmark_numbers = true;
        else
            % Delete them all
            for j=1:length(landmark_number_handles)
                delete(landmark_number_handles{j});
            end
            
            % Mark landmak numbers as not being shown
            showing_landmark_numbers = false;
        end
        
        
        
    end
    
    function flip_source_horizontal()
        %UNDOCUMENTED function that will flip the source image INCLUDING
        %THE IMAGE FILE ITSELF horizontally. Should actually be done during
        %preprocessing of the dataset, this function is for emergencies
        %only.
        
        source_image = fliplr(source_image);
        image_projection.CData = source_image;
        imwrite(source_image, paths.source);        
        
    end

    function auto_landmarks()
        %AUTO_LANDMARKS will attempt to automatically put landmarks on both
        %slices.
        
        
        
        
    end
end

